% DOES:
% * create non-movie time traces of different properties: mean fluorescence
% and number of nuclei transcribing grouped by eve stripe
% - detect locations of eve stripes / align different nuclei by stripe
% TODO:
% - allow averaging of multiple nuclei
% - calculate deviations from normal stripe positions
% - save .mat files to save repeat calculations?

%Dependencies: 
%loadCP.m, getParticlesInFrame.m, findStripes.m, plotByStripe.m,
%plotAllStripes.m, countStripes.m
%prefers export_fig.m to saveas.m

function [] = traceByStripe(varargin)
close all
set(0, 'defaulttextinterpreter', 'none'); %this isn't working idk 
%TODO - allow legend entries to have underscores w/o treating them as
%subscripts

%Control Flags
SaveFlag = 1;
AverageFlag = 0;
ClusterPlotFlag = ''; %This flag is used as an argument for findstripes

%Create folder to hold generated figures
if SaveFlag && ~exist([cd,filesep,'StripeFigures'], 'dir')
    mkdir([cd,filesep,'StripeFigures'])
end

%Parse inputs and determine the number of embryos we're working with
if ~isempty(varargin)
    %Loop through args to pull out numbers of embryos in each averaging
    %group.  n[] stores the number of each group to be averaged
    %length(n) is the number of different conditions being averaged
    %sum(n) is the total number of embryos being analyzed
    i = 1;
    while i <= length(varargin) && isnumeric(varargin{i})
        n(i) = varargin{i};
        i = i+1;
    end
    
    %keep track of total embryos / if there will be averaging
    nEmbryos = sum(n);
    nAvg = length(n);
    if nAvg > 1
        AverageFlag = 1;
    else
        clear
    end
    
    %parse options (loop through the rest of the args)
    for j = i:length(varargin)
        if strcmpi(varargin{j}, 'nosave')
            SaveFlag = 0;
        elseif strcmpi(varargin{j}, 'plotclusters')
            ClusterPlotFlag = 'Plot';
        end
    end
else
    %Default condition, 1 embryo, no averaging
    nEmbryos = 1;
end

display(['Running traceByStripe on ', num2str(nEmbryos), ' embryo(s)...'])

%Initialize data structures
E = cell(1,nEmbryos);
CP = cell(1,nEmbryos);
Prefix = cell(1, nEmbryos);
nParticles = NaN(1, nEmbryos);
nFrames = NaN(1, nEmbryos);
traceLength = 120; %number of frames
syncedFrames = NaN(nEmbryos, traceLength);
defaultStripes = 7;  %number of eve stripes (regions) we're looking at
stripes = cell(1,nEmbryos);

Embryos = struct; %main struct
Embryos(nEmbryos).stripeData = cell(traceLength, defaultStripes);
Embryos(nEmbryos).meanFluo = zeros(traceLength, defaultStripes);
Embryos(nEmbryos).activeNuclei = zeros(traceLength, defaultStripes);

%TODO - allow for infinite colors or cycle through or something
%Load relevant CompiledParticles data structures
for k = 1:nEmbryos
    loadPath = uigetdir(cd, ['Choose processed data set #', num2str(k)]);
    try
        E{k} = load([loadPath,filesep, 'CompiledParticles.mat'],...
            'CompiledParticles', 'ElapsedTime', 'nc14');
    catch
        error('CompiledParticles.mat not found')
    end
    
    %Get name of dataset
    dashes=strfind(loadPath,filesep);
    Prefix{k}=loadPath((dashes(end)+1):end);
    %Shorten name for convenience
    CP{k} = E{k}.CompiledParticles;
    %Get number of particles and length of each trace
    nParticles(k) = length(CP{k});
    nFrames(k) = length(E{k}.ElapsedTime);
    %Synchronize start point at nc14 mitosis, take [movieLength] frames
    syncedFrames(k,:) = (E{k}.nc14):(E{k}.nc14+traceLength-1);
    %TODO make number of frames customizable and/or fixed to length of
    %input traces.  also synchronize based on time and not frames
end

if nEmbryos > 1 && SaveFlag
    %Name to save figures
    filePrefix = input('Please enter a file prefix:\n?> ', 's');
    if AverageFlag && SaveFlag
        conditions = cell(1,nAvg);
        for j = 1:nAvg
            conditions{j} = input(...
                ['Name averaging condition ', num2str(j), '\n?>'], 's');
        end
    end
end

%--------------------------------------------------------------------------
%Begin main calculation loop: go through each embryo, through each time
%point, synced from the beginning of nc14 and through traceLength[] frames
%The syncedFrames[] array holds the corresponding frames for each trace

h = waitbar(0, 'Initializing');
for k = 1:nEmbryos
    %UI
    waitbar((k-1)/nEmbryos, h, ['Analyzing Embryo #', num2str(k),...
        ': Determining number of stripes in viewing window'])    
    fprintf('Analyzing Embryo %d/%d... ', k, nEmbryos)
    
    %Perform preanalysis - this will tell us how many stripes were captured
    %in the imaging window
    %This can be expanded to also determine a good transitionFrame and
    %startFrame for later
    nStripes = countStripes(CP{k}, syncedFrames(k,:));
    Embryos(k).stripeData = cell(traceLength, defaultStripes);

    fprintf('%d stripes detected\n', nStripes)
    fprintf('Clustering... ')
    %Run through time backwards to make initial stripe finding easier
    %Since the stripes are visible at the end of nc14 but not the
    %beginning, we can use clustering to identify stripes at later time
    %points but not earlier ones.  moving backward through time, we need to
    %change strategies at some point when clustering is no longer effective
    %i.e. the stripes are not distinguishable anymore, at which point we
    %will use the most recently calculated centroids to assign each nucleus
    %to a stripe in one step, without recalculating
    
    %Right now this is convoluted: start at startFrame, move back to zero,
    %then go forward from startFrame+1 until the end
    transitionFrame = 65;
    startFrame = 80;
    
    T = 0; %counter for loop/waitbar
    for t = [fliplr(1:startFrame), startFrame+1:traceLength]
        %Update waitbar
        T = T + 1;
        waitbar((traceLength*(k-1) + T)/(traceLength*nEmbryos), h,...
            ['Analyzing Embryo #', num2str(k),...
            ': Frame #', num2str(t), ' of nc14']);
        
        [pos, fluo] = getParticlesInFrame(CP{k}, syncedFrames(k,t));
        
        %STRIPE FINDING FLOW
        if t == startFrame || t == startFrame+1
            %First step, run findStripes with default centroids
            [Embryos(k).stripeData(t,:), centroids, stripes{k}] = ...
                findStripes(pos, fluo, nStripes, ClusterPlotFlag);
        elseif t > transitionFrame && t < startFrame
            %Frames with recognizable stripes, use clustering
            [Embryos(k).stripeData(t,:), centroids] = ...
                findStripes(pos, fluo, nStripes, centroids,...
                ClusterPlotFlag);
        else
            %Frames without recognizable stripes, use simple assignment
            Embryos(k).stripeData(t,:) = ...
                findStripes(pos, fluo, nStripes, centroids,...
                'AssignOnly', ClusterPlotFlag);
        end
                
        %Compute averages- still save the stripeData structure in case I
        %want to do something with it later, for now i'm only looking at
        %mean behavior.  TODO - standard deviations? stripe sizes?
        for s = 1:defaultStripes
            Embryos(k).meanFluo(t,s) = ...
                nanmean(Embryos(k).stripeData{t,s});
            Embryos(k).activeNuclei(t,s) = ...
                length(Embryos(k).stripeData{t,s});
        end
    end
    fprintf('Done!\n')
    
    %Plot and save
    [activeNucleiFig, meanFluoFig] = ...
        plotAllStripes(stripes{k}, Embryos(k), Prefix{k});
    if SaveFlag
        figSave(activeNucleiFig, [Prefix{k},'_activeNucleiFig'])
        figSave(meanFluoFig, [Prefix{k},'_meanFluoFig'])
    end
end
delete(h)

%--------------------------------------------------------------------------
%Plot each stripe on its own axis for multiple embryos, compare embryos
%stripe by stripe
if nEmbryos > 1;
    for s = 1:7
        stripeFig = plotByStripe(s,stripes,Embryos,Prefix);       
        %Save figure
        if SaveFlag
            figSave(stripeFig, [filePrefix, '_Stripe', num2str(s)]);
        end
    end
end

%Do averaging stuff here
if AverageFlag
    %create looping indices from n[] in terms of the embryos[k] so all the
    %embryos in condition j are numbered 
    %whichEmbryo(j) < k < whichEmbryo(j+1) - 1
    %e.g. if n = [2 5 3], whichEmbryo = [1 3 8 11];
    whichEmbryo = ones(1,nAvg+1);
    for j = 1:nAvg
        whichEmbryo(j+1:end) = ...
            whichEmbryo(j+1:end) + n(j)*ones(1,nAvg-j+1);
    end
    
    %hold new avg structures in a struct array, equivalent to Embryos()
    avgEmbryo = struct;
    %full stripe data will be necessary for errorbars, until then, w/e
    %avgEmbryo(nAvg).stripeData = cell(traceLength, defaultStripes);
    avgEmbryo(nAvg).meanFluo = zeros(traceLength, defaultStripes);
    avgEmbryo(nAvg).activeNuclei = zeros(traceLength, defaultStripes);
    
    %loop over averaging conditions
    for j = 1:nAvg       
        stripeSampleSize = NaN(1, defaultStripes);
        
        %loop over embryos in each averaging condition
        for k = whichEmbryo(j):whichEmbryo(j+1)-1
            %TODO normalize everything to have 7 stripes??
            %do averaging calculations
            avgEmbryo(j).meanFluo(:,stripes{k}) = ...
                avgEmbryo(j).meanFluo(:,stripes{k}) ...
                + Embryos(k).meanFluo;
            avgEmbryo(j).activeNuclei(:,stripes{k}) = ...
                avgEmbryo(j).activeNuclei(:,stripes{k}) ...
                + Embryos(k).activeNuclei;
            
            %how embryos many are we averaging over for each stripe?
            for s = 1:7
                if ~isempty(find(stripes{k}(stripes{k} == s),1))
                    stripeSampleSize(s) = stripeSampleSize(s) + 1;
                end
            end
        end
        
        %finish the averaging
        for s = 1:defaultStripes
            avgEmbryo(j).meanFluo(:,s) = ...
                avgEmbryo.meanFluo(:,s) / stripeSampleSize(s);
            avgEmbryo(j).activeNuclei(:,s) = ...
                avgEmbryo.activeNuclei(:,s) / stripeSampleSize(s);
        end
        %plot all stripes for this 'embryo'
        [activeNucleiFig, meanFluoFig] = ...
            plotAllStripes(1:7, avgEmbryo(j), conditions{j});
        %Save
        if SaveFlag
            figSave(activeNucleiFig, ...
                [conditions{j},'_activeNucleiFig'])
            figSave(meanFluoFig, ...
                [conditions{j},'_meanFluoFig'])
        end
    end
    
    %plot stripe by stripe comparisons
    for s = 1:7
        stripeFig = plotByStripe(s,1:7,avgEmbryo,Conditions);
        if SaveFlag
            figSave(stripeFig, [filePrefix, '_Stripe', s])
        end
    end
end

end

function [] = figSave(handle, name)
%Yes the export_fig syntax is weird but it works a lot better than
%MATLAB saveas, catch with saveas if export_fig not installed
try
    figure(handle)
    export_fig(sprintf('%s%sStripeFigures%s%s',cd,filesep,filesep,name));
catch
    warning(['Matlab saveas() sucks, try installing ',...
        'export_fig(): http://www.mathworks.com/matlabcentral',...
        '/fileexchange/23629-export-fig'])
    saveas(handle, ['StripeFigures', filesep, name]);
end
close handle
end