function [stripes, borders] = sortStripes(traces, nStripes)
%**************************************************************************
%
% Assign particles to stripes by drawing lines through the heatmap.
% Dependencies: none
%
%**************************************************************************

%Number of stripes exhibited by genotype - default is 7 but could be less
%in the case of stripe merging in mutants
if nargin < 2
    nStripes = 7;
end

stripes = [];

if isstruct(traces)
    beep;
    fprintf('Input must be an array of fluorescence time-traces')
    return
end

[traceLength, nParticles, ~] = size(traces);

%Bin traces to make heatmap for border specification
bins = 100;
binTraces = NaN(traceLength,bins);
[~,~, bin] = histcounts(nanmean(traces(:,:,3)),0:1/bins:1);

for t = 1:traceLength
    for x = 1:100
        %Total Fluorescence in bin
        binTraces(t,x) = nansum(traces(t,bin==x,1));
    end
end

%Allow user to specify each border on heatmap
heatmap = figure('Units', 'Normalized', 'Position', [0 0 1 1]);
imagesc(binTraces);
hold on;

borders = NaN(traceLength, nStripes+1);
for i = 1:nStripes+1
    figure(heatmap);
    title(['Draw border between stripes ',...
        num2str(i-1), ' and ', num2str(i)])
    [x, t] = ginputc('Color', [0.1 0.1 0.1], 'ShowPoints', true);
    %Check that selected points are legit, bounce back otherwise
    %**TODO
    
    %Interpolate drawn borders
    borders(:,i) = ...
        interp1([0; t; traceLength], [x(1); x; x(end)], 1:traceLength);
    
    %Plot line separating stripes
    plot(borders(:,i),1:traceLength, 'w--', 'LineWidth', 2)
end
close(heatmap)

%Assign particles to a stripe in each frame
stripes = NaN(traceLength,nParticles);
for t = 1:traceLength
    [~, ~, stripes(t,:)] = histcounts(traces(t,:,3), borders(t,:)/100);
end

%Examine traces that may have been split, resolve ambiguity of
%assignment (possibly split traces?)

%** TODO **

%Stripe assignment of a trace is the most frequent stripe assignment in
%each frame

stripes(stripes==0) = NaN;
stripes = mode(stripes);


%OLD VERSION BELOW, TAKES LONGER BUT IS MORE PRECISE. NOT QUITE FINISHED.
%**************************************************************************
%Use clustering to identify each stripe in an eve experiment, using the
%CompiledParticles structure.  Allow user to edit stripe assignments.
%-----> (Effectively a much better version of sortByStripe.m)
%Curation commands:
%
%   '1-9' (a number less than nStripes + 1)
%       - allows user to redifine the border between stripes n-1 and n
%       border between stripes 
%   '.' - advance one frame (last frame moves to first frame)
%   ',' - move back one frame (first frame moves to last frame)
%   'x' - exit & return sorted structure
%   'm' - shift all stripes left by one (i.e. stripes 2-7 become 1-6)
%   'n' - shift all strieps right by one
%   '+' - recluster current frame with one more stripe
%   '-' - recluster current frame with one less stripe
%
% Figures are only regenerated after valid inputs, to save frustrating 
% figure regeneration time
% %Filter out only nc14
% if isfield(cp, 'nc14')
%     cp.CompiledParticles = ...
%         cp.CompiledParticles([cp.CompiledParticles.FirstFrame] >= cp.nc14);
% else
%     warning('Field "nc14" does not exist, assuming all particles are nc14')
% end
% 
% %Get info about the length of traces (again only nc14)
% nFrames = length(cp.ElapsedTime(cp.nc14:end));
% frameShift = cp.nc14-1;         %Convert frame of movie to frame of nc14
% 
% %Figure out viewing window of the movie in terms of AP position
% DEFAULT_CENTROIDS = [0.32 0.40 0.46 0.525 0.59 0.69 0.78]';
% minAP = min([cp.CompiledParticles.APpos]);
% maxAP = max([cp.CompiledParticles.APpos]);
% 
% %Figure out which stripes we should be seeing based on window - this can be
% %edited by the use later but it gives a good starting point
% 
% %This is how we distinguish which stripes are actually visible
% whichStripes = DEFAULT_CENTROIDS > minAP & DEFAULT_CENTROIDS < maxAP;
% 
% %Use clustering to specify stripes initially 
% centroids = NaN(nFrames, nStripes);
% clusterRange = 46:85;
% earlyRange = 1:clusterRange(1)-1;
% lateRange = clusterRange(end)+1:nFrames;
% 
% APpos = cell(nFrames, 1);   %Save all AP- and y- positions in parallel 
% yPos = cell(nFrames, 1);    %cells to save time during user interaction
% for t = 1:nFrames
%     [~,~, APpos{t}, yPos{t}] = particlesInFrame(...
%         cp.CompiledParticles, t + frameShift, 'MedianAP', 'yPos');
%     
%     %Do clustering in the mid-nc14 period while stripes are well-defined
%     if any(clusterRange == t)
%         [~, centroids(t,whichStripes)] = ...
%             kmeans(APpos{t}', sum(whichStripes),...
%             'Start', DEFAULT_CENTROIDS(whichStripes));  
%         centroids(t,whichStripes) = sort(centroids(t,whichStripes));
%     end
% end
% 
% %Set the centroids for the edges
% earlyCentroids = mean(centroids(clusterRange(1:3),:));
% centroids(earlyRange,:) = repmat(earlyCentroids,length(earlyRange),1);
% 
% lateCentroids = mean(centroids(clusterRange(end-2:end),:));
% centroids(lateRange,:) = repmat(lateCentroids,length(lateRange),1);
% 
% %Initialize borders
% borders = NaN(nFrames, nStripes+1);
% for t = 1:nFrames
%     %Take midpoints of centroids as borders for now + edges
%     borders(t,:) = [0.2 conv(centroids(t,:),[0.5 0.5],'valid') 0.85];
%     if isnan(borders(t,2))
%         borders(t,1) = NaN;
%         borders(t,2) = 0.2;
%     end
%     
%     if isnan(borders(t,end-1))
%         borders(t,end) = NaN;
%         borders(t,end-1) = 0.85;
%     end
% end

%---------------------MANUAL CURATION STARTS HERE--------------------------

% %Counter to keep track of where we are
% t = 45;
% 
% %Initialize figures
% bg = [0.1 0.1 0.1]; %Color for figure background
% fg = [0.9 0.9 0.9]; %Color for figure foreground
% particleFig = ...
%     figure('Units', 'Normalized', 'Position', [0 0.5 1, 0.5],...
%     'Color', bg, 'MenuBar', 'None', 'ToolBar', 'None',...
%     'NumberTitle', 'off');%, 'WindowStyle', 'Modal');
% histFig = ...
%     figure('Units', 'Normalized', 'Position', [0 0 0.6, 0.4],...
%     'MenuBar', 'None', 'ToolBar', 'None', 'NumberTitle', 'off');
% % heatFig = ...
% %     figure('Units', 'Normalized', 'Position', [0.6 0 0.4 0.4],...
% %     'MenuBar', 'None', 'ToolBar', 'None', 'NumberTitle', 'off');

% while true
%     % **** Generate particleFig ****
%     
%     %As a display tool, we use APpos as the x position (so that we don't
%     %have to worry about any conversion) and yPos as y Position because it
%     %isn't relevant for defining the borders of stripe expression
%     
%     %Display each particle in a different color based on stripe assignment
%     figure(particleFig)
%     cla; hold on;
%     for s = 1:nStripes
%         particlesInStripe = APpos{t} > borders(t,s) ...
%             & APpos{t} < borders(t,s+1);
%         plot(APpos{t}(particlesInStripe), yPos{t}(particlesInStripe), 'o');
%         text(centroids(t,s), 245, num2str(s),...
%             'Color', fg, 'FontWeight', 'bold', 'EdgeColor', fg)
%     end
%     
%     %Display lines separating each stripe
%     lines = cell(1,nStripes+1);
%     for b = 1:nStripes+1
%         lines{b} = plot(borders(t, b)*ones(1,2), [0 256],...
%             '--', 'Color', fg);
%     end
% 
%     %Specify axis properties etc.
%     axis([0.2 0.85 0 256])
%     particleFig.Name = ['Frame ', num2str(t), ' of nc14'];
%     ax = gca;
%     ax.Color = bg; ax.XColor = fg; ax.YColor = 'none';
%     
%     % **** generate histFig ****
%     
%     figure(histFig)
%     cla; hold on;
%     for s = 1:nStripes
%         if ~isnan(borders(t,s) + borders(t,s+1))
%             %Histogram of particles per bin (per stripe)
%             [posCounts, edges] = histcounts(...
%                 APpos{t}, borders(t,s):0.01:borders(t,s+1));
%             centers = [edges(1)-0.005, edges + 0.005];
%             plot(centers, [0, posCounts, 0], 'Linewidth', 2)
%             
%             %Centroid of each stripe (Repeat same color)
%             ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
%             plot([centroids(t, s) centroids(t, s)], [0 100], '--')
%         end
%     end
%     axis([0.2 0.85 0 10])
%     
%     % **** generate heatFig ****
% %     
% %     figure(heatFig)
% %     cla; hold on;
%     % **** User control ****
%     
%     %Get user input; do nothing unless input key is valid (see above)
%     %Set dummy value for input key to enter loop that waits for input
%     key = [];
%     while isempty(key) || any(key == 'ij') || ...
%             (~(str2double(key) <= nStripes+1) && all(key ~= 'x.,mn-+'))
%         figure(particleFig)
%         waitforbuttonpress;
%         key = particleFig.CurrentCharacter;
%     end
%     
%     %Exit program and save results
%     if key == 'x'
%         break
%         
%     %Change frame, use wraparound
%     elseif key == '.'
%         if t < traceLength
%             t = t + 1;
%         else
%             t = 1;
%         end
%     elseif key == ','
%         if t > 1
%             t = t - 1;
%         else
%             t = traceLength;
%         end
%         
%     %Shift stripes left or right
%     elseif key == 'm'
%         if isnan(centroids(t,end))
%             centroids(t,2:end) = centroids(t,1:end-1);
%             borders(t,2:end) = borders(t,1:end-1);
%             centroids(t,1) = NaN;
%             borders(t,1) = NaN;
%         else
%             beep;
%             warning('Cannot shift right')
%         end
%     elseif key == 'n'
%         if isnan(centroids(t,1))
%             centroids(t,1:end-1) = centroids(t,2:end);
%             borders(t,1:end-1) = borders(t,2:end);
%             centroids(t,end) = NaN;
%             borders(t,end) = NaN;
%         else
%             beep;
%             warning('Cannot shift left')
%         end
%         
%     %Add stripe
%     elseif key == '+'
%         nStripesInFrame = sum(~isnan(centroids(t,:)));
%         if  nStripesInFrame < nStripes
%             %Recluster
%             [~, newCentroids] = kmeans(APpos{t}', nStripesInFrame+1);
%             newCentroids = sort(newCentroids);
%             centroids(t,:) = NaN(1,nStripes);
%             centroids(t,1:nStripesInFrame+1) = newCentroids;
%             
%             %Redefine borders
%             borders(t,:) = ...
%                 [0.2 conv(centroids(t,:),[0.5 0.5],'valid') 0.85];
%         end
%         
%     %Border redefinition
%     else
%         %Input is a number, parse and do stuff with it
%         thisBorder = str2double(key);
%         lines{thisBorder}.Visible = 'off';
%         
%         %Allow for user to input new stripe location
%         oldPos = borders(t,thisBorder);
%         [newBorder, ~] = ginput(1);
%         %[newBorder, ~] = ginputc(1,...
%         %    'Color', fg, 'FigHandle', particleFig);
%         
%         %Catch attempted reassignments beyond other borders
%         if (thisBorder ~= nStripes+1 && ...
%                 newBorder > borders(t, thisBorder+1)) || ...
%                 (thisBorder ~= 1 && newBorder < borders(t, thisBorder-1))
%             beep;
%             warning('Cannot move border beyond adjacent borders')
%             continue
%         end
%         
%         borders(t,thisBorder) = newBorder;
%         
%         %Output results to command line
%         fprintf('Border %i in frame %i redefined.\n', thisBorder, t);
%         fprintf('Old position: %f  New Position: %f\n\n',...
%             oldPos, borders(t,thisBorder))
%         
%         %Recalculate centroids
%         for s = 1:nStripes
%             particlesInStripe = APpos{t} > borders(t,s) ...
%                 & APpos{t} < borders(t,s+1);
%             centroids(t,s) = mean(APpos{t}(particlesInStripe));
%         end
%     end
% end
% close('all', 'force')


end