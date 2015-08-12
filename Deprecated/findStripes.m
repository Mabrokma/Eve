%use 1-D kmeans clustering for finding stripes

function [stripeData, nucPerStripe, centroids, stripes] = ...
    findStripes(APpos, fluo, stripes, nNuclei, varargin)
%Flags
AssignOnlyFlag = 0;
PlotFlag = 0;

DEFAULT_CENTROIDS = [0.3 0.38 0.46 0.525 0.59 0.67 0.76]';
MAX_STRIPE_RADIUS = 0.04;

if length(stripes) == 1;
    nStripes = stripes;
else
    nStripes = length(stripes);
end

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
    %The window size can't always track all stripes for any given embryo
    %so we need to choose what to leave out, do this by looking at
    %max and min of APpos and tossing the farther from the edge bin
    
    if nStripes == 6;
        skipFirst = abs(min(APpos) - DEFAULT_CENTROIDS(1)) > ...
            abs(max(APpos) - DEFAULT_CENTROIDS(7));
        if skipFirst
            centroids = DEFAULT_CENTROIDS(2:7);
            stripes = 2:7;
        else
            centroids = DEFAULT_CENTROIDS(1:6);
            stripes = 1:6;
        end
    elseif nStripes == 7;
        centroids = DEFAULT_CENTROIDS;
        stripes = 1:7;
    elseif nStripes == 5;
        centroids = DEFAULT_CENTROIDS(2:6);
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
        kmeans(APpos', nStripes, 'Start', centroids, 'MaxIter', 20);
end

%filter outliers (d > MAX_STRIPE_RADIUS)
for i = 1:length(APpos)
    if APpos(i) - centroids(whichCluster(i)) > MAX_STRIPE_RADIUS
        whichCluster(i) = NaN;
    end
end


%convert everything to terms of absolute stripes---------------------------
%First sort clusters in order and assign whichStripe as a sorted
%whichCluster
[centroids, stripeToCluster] = sort(centroids);

whichStripe = NaN(1,length(whichCluster));
for s = 1:length(centroids)
    whichStripe(whichCluster==stripeToCluster(s)) = s;
end

nucWhichStripe = NaN(1,length(nNuclei));
for i = 1:length(nNuclei)
    bestD = Inf;
    %cycle through centroids, find out which one is closest
    for j = 1:nStripes
        if abs(centroids(j) - (i-1)/100) < bestD;
            bestD = abs(centroids(j) - (i-1)/100);
            if bestD < MAX_STRIPE_RADIUS;
                nucWhichStripe(i) = j;
            end
        end
    end
end


%If stripe 1 is not the first stripe, shift everything by 1
%todo - deal with fewer stripes than 5
if stripes(1) == 2
    whichStripe = whichStripe + 1;
    nucWhichStripe = nucWhichStripe + 1;
end

%Calculate total nuclei in each stripe
nucPerStripe = zeros(1,7);
for i = 1:length(nNuclei)
    try
        nucPerStripe(nucWhichStripe(i)) = ...
            nucPerStripe(nucWhichStripe(i)) + nNuclei(i);
    catch
        %Do nothing if nucWhichStripe(i) == NaN
    end
end
%Convert zeros to nans
nucPerStripe(nucPerStripe == 0) = NaN;

%Assign fluorescence values to the relevant stripe
stripeData = cell(1,7);
APbyStripe = cell(1,7);
for i = 1:length(APpos)
    try
        stripeData{whichStripe(i)}(end+1) = fluo(i);
        APbyStripe{whichStripe(i)}(end+1) = APpos(i);
    catch
        %This is for entries in whichStripe that were replaced with NaNs
    end
end

%Plot results for debugging------------------------------------------------
if PlotFlag
    ClusterFig = figure;
    
    hold on
    for s = 1:nStripes
        [posCounts, edges] = histcounts(APbyStripe{s}, 'BinWidth', 0.01);
        centers = [edges(1)-0.005, edges + 0.005];
        plot(centers, [0, posCounts, 0], 'Linewidth', 2)
        %Repeat same color
        ax = gca; ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
        
        plot([centroids(s) centroids(s)], [0 100],...
            '--', 'Linewidth', 2)
    end
    axis([0.2 0.85 0 10])
    
    pause(0.25)
    close(ClusterFig)
end