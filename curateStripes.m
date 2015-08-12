function [centroids, borders] = curateStripes(cp, centroids)
%**************************************************************************
%Curate the borders between eve stripes in a manner similar to
%CheckParticleTracking()
%Code responds to:
%
%   '(n)' (a number less than nStripes + 1)
%       - allows user to redifine the border between stripes n-1 and n
%       border between stripes 
%   '.' - advance one frame (last frame moves to first frame)
%   ',' - move back one frame (first frame moves to last frame)
%   'x' - exit
%
%And otherwise does nothing.  Figures are only regenerated after valid
%inputs, to save frustrating figure regeneration time
%**************************************************************************


%Initialize borders
[traceLength, nStripes] = size(centroids);
borders = NaN(traceLength, nStripes+1);
for t = 1:traceLength
    %Take midpoints of centroids as borders for now + edges
    borders(t,:) = [0.2 conv(centroids(t,:),[0.5 0.5],'valid') 0.85];
end

%Counter to keep track of where we are
frameShift = cp.nc14 - 1;
currentFrame = 45;

%Initialize figures
bg = [0.1 0.1 0.1]; %Color for figure background
fg = [0.9 0.9 0.9]; %Color for figure foreground
particleFig = ...
    figure('Units', 'Normalized', 'Position', [0 0.5 1, 0.5],...
    'Color', bg, 'MenuBar', 'None', 'ToolBar', 'None',...
    'NumberTitle', 'off', 'WindowStyle', 'Modal');%,...
    %'CloseRequestFcn', 'return;');
histFig = ...
    figure('Units', 'Normalized', 'Position', [0 0 1, 0.45],...
    'MenuBar', 'None', 'ToolBar', 'None', 'NumberTitle', 'off');

while true
    %Get particles in frame, assign to stripes
    %As a display tool, we use APpos as the x position (so that we don't
    %have to worry about any conversion) and yPos as y Position because it
    %isn't relevant for defining the borders of stripe expression
    [~,~, APpos, yPos] = particlesInFrame(...
        cp.CompiledParticles, currentFrame + frameShift, 'APpos', 'yPos');
    
%generate particleFig------------------------------------------------------
    %Display each particle in a different color based on stripe assignment
    figure(particleFig)
    cla; hold on;
    for s = 1:nStripes
        particlesInStripe = APpos > borders(currentFrame,s) ...
            & APpos < borders(currentFrame,s+1);
        plot(APpos(particlesInStripe), yPos(particlesInStripe), 'o');
        text(centroids(currentFrame,s), 245, num2str(s),...
            'Color', fg, 'FontWeight', 'bold', 'EdgeColor', fg)
    end
    
    %Display lines separating each stripe
    lines = cell(1,nStripes+1);
    for b = 1:nStripes+1
        lines{b} = plot(borders(currentFrame, b)*ones(1,2), [0 256],...
            '--', 'Color', fg);
    end

    %Specify axis etc.
    axis([0.2 0.85 0 256])
    particleFig.Name = ['Frame ', num2str(currentFrame), ' of nc14'];
    ax = gca;
    ax.Color = bg;
    ax.YColor = 'none';
    ax.XColor = fg;
    
%generate histFig----------------------------------------------------------
    figure(histFig)
    cla; hold on;
    for s = 1:nStripes
        %Histogram of particles per bin (per stripe)
        [posCounts, edges] = histcounts(...
            APpos, borders(currentFrame,s):0.01:borders(currentFrame,s+1));
        centers = [edges(1)-0.005, edges + 0.005];
        plot(centers, [0, posCounts, 0], 'Linewidth', 2)
        %Centroid of each stripe (Repeat same color)
        ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
        plot([centroids(currentFrame, s) centroids(currentFrame, s)],...
            [0 100],'--')
    end
    axis([0.2 0.85 0 10])
    
%Control-------------------------------------------------------------------
    %Get user input; do nothing unless input key is valid (see above)
    %Set dummy value for input key to enter loop that waits for input
    key = ' ';
    while ~(str2double(key) <= nStripes+1) ...
            && key ~= 'x' && key ~= '.' && key ~= ',' ...
            || key == 'i' || key == 'j'
        figure(particleFig)
        waitforbuttonpress;
        key = particleFig.CurrentCharacter;
    end
    
    %Exit program and save results
    if key == 'x'
        break
        
    %Change frame, use wraparound
    elseif key == '.'
        if currentFrame < traceLength
            currentFrame = currentFrame + 1;
        else
            currentFrame = 1;
        end
    elseif key == ','
        if currentFrame > 1
            currentFrame = currentFrame - 1;
        else
            currentFrame = traceLength;
        end
    
    %Actually do border redefinition
    else
        %Input is a number, parse and do stuff with it
        thisBorder = str2double(key);
        lines{thisBorder}.Visible = 'off';
        
        %Allow for user to input new stripe location
        oldPos = borders(currentFrame,thisBorder);
        [borders(currentFrame,thisBorder), ~] = ginput(1);
        
        %Output results to command line
        fprintf('Border %i in frame %i redefined.\n',...
            thisBorder, currentFrame);
        fprintf('Old position: %d  New Position: %d\n\n',...
            oldPos, borders(currentFrame,thisBorder))
        
        %Recalculate centroids
        % * * TODO
    end
end

close all
end


function countStripes()
potentialCentroids = cell(1, length(clusterRange));
foundStripes = NaN(length(clusterRange),1);
left = NaN(length(clusterRange),1);
APpos = cell(1,length(clusterRange));
binSize = cp.APbinID(2);

%NOTE: APpos is indexed from 1:length(clusterRange) but the values of
%clusterRange do not begin at 1; t corresponds to frame, T to index
T = 0;
for t = clusterRange
    %Use capital T for indexing 
    T = T+1;
    APpos{T} = getParticlesInFrame(cp.CompiledParticles, frameOfnc14(t));
    
    %Make histogram of frame, smooth, count peaks as proxy for stripes
    posCounts = histcounts(APpos{T}, cp.APbinID);
    w = 4; %moving window size
    smoothCounts = conv(posCounts, ones(1,w)/w, 'same'); %filter
    [~,potentialCentroids{T}] = findpeaks(smoothCounts);
    potentialCentroids{T} = potentialCentroids{T}*binSize - binSize/2;
    foundStripes(t) = length(potentialCentroids{T});
    
    distancesFromCentroids = ...
        abs(potentialCentroids{T}(1)-DEFAULT_CENTROIDS);
    %Compare found centroids to default centroids
    [d, left(T)] = min(distancesFromCentroids);
end
nStripes = mode(foundStripes);
firstStripe = mode(left);
lastStripe = firstStripe+nStripes-1;
stripeShift = firstStripe - 1;
end

function assignToStripe()

end

function calculateCentroids(borders)

end