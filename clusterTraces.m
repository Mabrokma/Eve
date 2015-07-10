%Cluster individual traces and determine the geographic proximity of
%particles with similar flourescence traces

function [results] = clusterTraces(trace, meanAP, nClusters)
%How should NaNs be handled?  Method 1: replace with mean-3std
trace(isnan(trace)) = 0;

whichCluster = kmeans(trace', nClusters);
[~, order] = sort(whichCluster);

imagesc(trace(:,order))

results = zeros(nClusters,3);
for i = 1:nClusters
    results(i,1) = i;
    results(i,2) = mean(meanAP(whichCluster==i));
    results(i,3) = std(meanAP(whichCluster==i));
end
