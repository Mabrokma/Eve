%use 1-D kmeans clustering for finding stripes

function [stripeData, centroids, stripes] = ...
    findStripes(APpos, fluo, nStripes, varargin)
%Flags
AssignOnlyFlag = 0;
PlotFlag = 0;

if ~isempty(varargin)
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'assignonly')
            AssignOnlyFlag = 1;
        elseif strcmpi(varargin{i}, 'plot')
            PlotFlag = 1;
        end
    end
end

%See if centroid seeds were provided as an argument
if ~isempty(varargin) && isnumeric(varargin{1})
    centroids = varargin{1};
    nStripes = length(centroids);
else
    defaultCentroids = [0.3 0.38 0.46 0.525 0.59 0.67 0.76]';
    
    %The window size can't always track all stripes for any given embryo
    %so we need to choose what to leave out, do this by looking at
    %max and min of APpos and tossing the farther from the edge bin
    
    if nStripes == 6;
        skipFirst = abs(min(APpos) - defaultCentroids(1)) > ...
            abs(max(APpos) - defaultCentroids(7));
        if skipFirst
            centroids = defaultCentroids(2:7);
            stripes = 2:7;
        else
            centroids = defaultCentroids(1:6);
            stripes = 1:6;
        end
    elseif nStripes == 7;
        centroids = defaultCentroids;
        stripes = 1:7;
    elseif nStripes == 5;
        centroids = defaultCentroids(2:6);
        stripes = 2:6;
    else
        error('Stripes too irregular to cluster w/o seeds')
    end
    
    %maybe initialize centroids some other way besides defaults??!?!?
    %idea: find maxima of filtered histograms; potential filter schemes
    % = moving window average, spline?
    
    %this is sorta done by countStripes() now but there is no communication
    %between the two functions
end

%Determine mode and do clustering/assignment-------------------------------
if AssignOnlyFlag
    %Simple assignment mode, just assign particles to nearest centroid
    whichCluster = zeros(length(APpos), 1);
    for i = 1:length(APpos)
        bestD = Inf;
        %cycle through centroids, find out which one is closest
        for j = 1:nStripes
            if abs(centroids(j) - APpos(i)) < bestD;
                bestD = abs(centroids(j) - APpos(i));
                whichCluster(i) = j;
            end
        end
    end
else
    %Clustering mode: Use built-in matlab kmeans to cluster particles
    %(inputs need to be column vectors)
    [whichCluster, centroids] = ...
        kmeans(APpos', nStripes, 'Start', centroids, 'MaxIter', 10);
end

%convert everything to terms of absolute stripes
stripeToCluster = zeros(1,nStripes);
whichStripe = zeros(1,length(whichCluster));
for i = 1:nStripes
    bestD = inf;
    %compare centroids to the default centroids, assign based on the
    %default stripe centroid to which it's closest
    for s = 1:7
        thisD = abs(centroids(i) - defaultCentroids(s));
        if thisD < bestD
            bestD = thisD;
            stripeToCluster(s) = i;
        end
    end
    %take whichcluster and turn it into whichstripe
    %find locations of each centroid assignment
    whichStripe(find(whichCluster==stripeToCluster(s))) = s;
end

%Assign fluorescence values to the relevant stripe
stripeData = cell(1,7);
APbyStripe = cell(1,7);
for i = 1:length(APpos)
    stripeData{whichStripe(i)}(end+1) = fluo(i);
    APbyStripe{whichStripe(i)}(end+1) = APpos(i);
end

%Plot results for debugging------------------------------------------------
if PlotFlag
    ClusterFig = figure;
    
    hold on
    for i = 1:nStripes
        [posCounts, edges] = histcounts(APbyStripe{i}, 'BinWidth', 0.01);
        centers = [edges(1)-0.005, edges + 0.005];
        plot(centers, [0, posCounts, 0], 'Linewidth', 2)
        plot(centroids(i)*ones(1,2), [0 100],'--', 'Linewidth', 2)
    end
    axis([0.25 0.75 0 10])
    
    pause(0.25)
    close(ClusterFig)
end