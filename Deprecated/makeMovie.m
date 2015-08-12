%%Export time traces of certain particle properties to movies
% RW begun 6/19/15

%Large scale goals / TODO:
% - allow for averaging of multiple embryos
% - detect / align eve stripes
% - allow specification of nuclear cycle
% - make movies of arbitrary data in the CompiledParticles structure
% - add error handling / prevent crashing


function [] = makeMovie(varargin)
close all
%Parse inputs, determine averaging vs not averaging
l = length(varargin);
if l == 0
    nEmbryos = 1;
elseif l == 1 && isnumeric(varargin{1})
    nEmbryos = varargin{1};
elseif l == 1 && ischar(varargin{1})
    error('Not done with manual folder specification')
elseif l == 2 && isnumeric(varargin{2})
    error('Not done with averaging yet')
else
    error('You goofed, fix the args')
end

%Initialize data structures
E = cell(1,nEmbryos);
C = cell(1,nEmbryos);
Prefix = cell(1, nEmbryos);
nParticles = NaN(1, nEmbryos);
nFrames = NaN(1, nEmbryos);
movieLength = 120; %number of frames
syncedFrames = NaN(nEmbryos, movieLength);
%TODO - allow for infinite colors or cycle through or something

%Properties to use for distributions (binning etc.)
%TODO - make this customizable via input
%TODO - make axes depend on these
posBins = 0.25:0.01:0.75;
posBinCenters = posBins + 0.005; 
posBinCenters(end) = [];

fluoBins = 0:250:7000;
fluoBinCenters = fluoBins + 125; 
fluoBinCenters(end) = [];

%Get relevant directories
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
    C{k} = E{k}.CompiledParticles;
    
    %Get number of particles and length of each trace
    nParticles(k) = length(C{k});
    nFrames(k) = length(E{k}.ElapsedTime);
    
    %Synchronize start point as beginning of nc14, take [movieLength] frames
    syncedFrames(k,:) = E{k}.nc14:E{k}.nc14+movieLength-1;
    %TODO make number of frames customizable and/or fixed to length of
    %input traces.  also synchronize based on time and not frames
end

%initialize VideoWriter objects
if nEmbryos > 1
    filePrefix = input('Please enter a file prefix:\n?> ', 's');
else
    filePrefix = Prefix{1};
end

fluoDistMovie = VideoWriter([filePrefix, '_fluoDistMovie']);
fluoDistMovie.FrameRate = 5;
open(fluoDistMovie);

posDistMovie = VideoWriter([filePrefix, '_posDistMovie']);
posDistMovie.FrameRate = 5;
open(posDistMovie);

%open figure windows
fluoDistFig = figure;
posDistFig = figure;

for t = 1:movieLength
    cla; %Clear previous plot
    hold on
    
    %Reset histogram counting structures
    allPosInFrame = cell(1,nEmbryos);
    allFluoInFrame = cell(1,nEmbryos);
    %Run searches on each embryo
    %TODO - allow for averaging of multiple embryos
    %TODO - fix distributions to represent fractions of (approved) nuclei
    % i.e. transcribing nuclei as a fraction of total nuclei in a given
    % slice, excluding disapproved nuclei (?) - see if this makes a dif
    for k = 1:nEmbryos
        [allPosInFrame{k}, allFluoInFrame{k}] = ...
            getParticlesInFrame(C{k}, syncedFrames(k,t));
    end
%--------------------------------------------------------------------------

    %Make fluo histogram
    figure(fluoDistFig);
    cla; %Clear previous plot
    hold on
    
    fluoCounts = cell(1, nEmbryos);
    h = zeros(1,nEmbryos);
    for k = 1:nEmbryos
        fluoCounts{k} = histcounts(allFluoInFrame{k}, fluoBins);
        meanFluo = mean(allFluoInFrame{k});
        h(k) = plot(fluoBinCenters, fluoCounts{k});
        plot([meanFluo meanFluo], [0 25],...
            '--', 'LineWidth', 2)
    end
    
    %Customize axes for fluoDistFig
    axis([0, 7000, 0, 25])
    title(['Frame of nc14: ', num2str(t)])
    xlabel('Fluorescence Intensity')
    ylabel('Frequency')
    for k = 1:nEmbryos
        legend(h, Prefix, 'Location', 'NorthEast')
    end

    %Set area of plot to capture
    coords = fluoDistFig.Position;
    rect = [0, 0, coords(3), coords(4)]; 

    %Capture frame and write to movie
    movieFrame = getframe(fluoDistFig, rect);
    writeVideo(fluoDistMovie, movieFrame);
%--------------------------------------------------------------------------
    %Make pos histogram
    figure(posDistFig);
    cla; %Clear previous plot
    hold on
    
    posCounts = cell(1,nEmbryos);
    for k = 1:nEmbryos
        posCounts{k} = histcounts(allPosInFrame{k}, posBins);
        plot(posBinCenters, posCounts{k},...
            'Linewidth', 2)
    end
    
    %Customize axes for posDistFig
    axis([0.25, 0.75, 0, 15])
    title(['Frame of nc14: ', num2str(t)])
    xlabel('AP Position')
    ylabel('Frequency')
    legend(Prefix, 'Location', 'NorthEast')
    %Set area of plot to capture
    coords = posDistFig.Position;
    rect = [0, 0, coords(3), coords(4)]; 

    %Capture frame and write to movie
    movieFrame = getframe(posDistFig, rect);
    writeVideo(posDistMovie, movieFrame);
end

%Close figures
close all
%Close movies
close(fluoDistMovie);
close(posDistMovie);