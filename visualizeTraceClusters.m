function visualizeTraceClusters(rawTraces, clusters, time14)
%**************************************************************************
%Allow for user to manually click through a virtual movie with the
%particles colored according to clusters[].  clusters[] should have the
%same length as the number of columns of rawTraces, and time14 is the
%portion of ElapsedTime that corresponds to nc14.
%
% In the process of allowing one to cycle through different numbers of
% clusters, in which case clusters[] would be a cell array or 2D vector
% containing cluster assignments for each k
%
%Controls:
%
%   '.' - advance one frame
%   ',' - reverse one frame
%   'm' - increase number of clusters
%   'n' - decrease number of clusters
%   'x' - exit animation
%   's' - generate/save movie (loop through all t and use videowriter
%   object to save each figure)
%
% Dependencies: none
% RW 8/2015
%**************************************************************************

%Extract relevant parameters from raw traces
nFrames = size(rawTraces, 1);
nClusters = max(clusters);

%Make other plots showing behavior of clusters over time
clusterFig = figure('Units', 'Normalized', 'Position', [0 0 0.5 0.4]);
[~, order] = sort(clusters);
%Plot each trace by cluster assignment
imagesc(rawTraces(:,order,1)); colormap jet;
title(['Individual Traces: ' num2str(nClusters), ' clusters'])
%Add lines separating clusters and labels
for i = 1:nClusters
        hold on
        %Label each cluster
        text(mean(find(sort(clusters)==i)), 5, num2str(i),...
            'Color', [1 1 1], 'FontWeight', 'bold', 'EdgeColor', [1 1 1]);
        
        %Plot lines between clusters
        if i ~= nClusters
            plot(find(sort(clusters)==i, 1, 'last')*ones(1,2),...
                [0 100], 'm', 'Linewidth', 2)
        end
        ax = gca;
        ax.XTick = [];
        ax.XTickLabel = [];
        
        xlabel('Cluster')
        ylabel('Time')
end

%Plot mean traces for each cluster
traceFig = figure('Units', 'Normalized', 'Position', [0.5 0 0.5 0.4]);
leg = cell(1,nClusters);
fluoTrace = rawTraces(:,:,1);
fluoTrace(isnan(fluoTrace)) = 0;
for i = 1:nClusters
    hold on
    plot(time14, mean(fluoTrace(:,clusters == i),2),...
        'LineWidth', 2);
    xlabel('Frame Number')
    ylabel('Mean Fluorescence')
    leg{i} = ['Cluster ', num2str(i)];
end
legend(leg, 'Location', 'Best')


%Enter user controlled visualization of individual particles 
t = 1;
particleFig = plotFrame(t, clusters, rawTraces);
while true
    %Get user input; do nothing unless input key is valid (see above)
    %Set dummy value for input key to enter loop that waits for input
    key = [];
    while isempty(key) || all(key ~= '.,xs')
        figure(particleFig)
        waitforbuttonpress;
        key = particleFig.CurrentCharacter;
    end
    
    %Exit program and save results
    if key == 'x'
        break
        
    %Change frame, use wraparound
    elseif key == '.'
        if t < nFrames
            t = t + 1;
        else
            t = 1;
        end
    elseif key == ','
        if t > 1
            t = t - 1;
        else
            t = nFrames;
        end
    elseif key == 's'
        warning('Workinonit')
        break
    end
    
    particleFig = plotFrame(t, clusters, rawTraces, particleFig);
end

close all

end


function particleFig = plotFrame(t,clusters, rawTraces, particleFig)
bg = [0.1 0.1 0.1]; %Color for figure background
fg = [0.9 0.9 0.9]; %Color for figure foreground
if nargin < 4
    %Initialize figure
    particleFig = ...
        figure('Units', 'Normalized', 'Position', [0 0.55 1, 0.45],...
        'Color', bg, 'MenuBar', 'None', 'ToolBar', 'None',...
        'NumberTitle', 'off');%, 'WindowStyle', 'Modal');
else
    figure(particleFig)
end

nClusters = max(clusters);

%Display each particle in a different color based on cluster assignment
cla; hold on;
whichParticles = cell(1,nClusters);
for i = 1:nClusters
    whichParticles{i} = (clusters == i);
    plot(rawTraces(t, whichParticles{i}, 3),...
        rawTraces(t, whichParticles{i}, 5),'o');
end

%Specify axis properties etc.
axis([0.2 0.85 0 256])
title(['Active nuclei colored by cluster; Frame: ', num2str(t)])
particleFig.Name = ['Frame ', num2str(t), ' of nc14'];
ax = gca;
ax.Color = bg; ax.XColor = fg; ax.YColor = 'none';

particleFig = gcf;
end

function [clusterFig traceFig] = ...
    plotCluster(k, clusters, rawTraces, time14)
%Make other plots showing behavior of clusters over time
clusterFig = figure('Units', 'Normalized', 'Position', [0 0 0.5 0.4]);
[~, order] = sort(clusters);
%Plot each trace by cluster assignment
imagesc(rawTraces(:,order,1)); colormap jet;
title(['Individual Traces: ' num2str(nClusters), ' clusters'])
%Add lines separating clusters and labels
for i = 1:nClusters
        hold on
        %Label each cluster
        text(mean(find(sort(clusters)==i)), 5, num2str(i),...
            'Color', [1 1 1], 'FontWeight', 'bold', 'EdgeColor', [1 1 1]);
        
        %Plot lines between clusters
        if i ~= nClusters
            plot(find(sort(clusters)==i, 1, 'last')*ones(1,2),...
                [0 100], 'm', 'Linewidth', 2)
        end
        ax = gca;
        ax.XTick = [];
        ax.XTickLabel = [];
        
        xlabel('Cluster')
        ylabel('Time')
end

%Plot mean traces for each cluster
traceFig = figure('Units', 'Normalized', 'Position', [0.5 0 0.5 0.4]);
leg = cell(1,nClusters);
fluoTrace = rawTraces(:,:,1);
fluoTrace(isnan(fluoTrace)) = 0;
for i = 1:nClusters
    hold on
    plot(time14, mean(fluoTrace(:,clusters == i),2),...
        'LineWidth', 2);
    xlabel('Frame Number')
    ylabel('Mean Fluorescence')
    leg{i} = ['Cluster ', num2str(i)];
end
legend(leg, 'Location', 'Best')
end