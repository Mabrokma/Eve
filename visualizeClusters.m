function visualizeClusters(clusters, rawTraces, time)
%**************************************************************************
%Allow for user to manually click through a virtual movie with the
%particles colored according to clusters{}.  Each matrix in clusters{} 
%should have the same length as the number of columns of rawTraces.
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
% Dependencies: clusterTraces
% RW 8/2015
%**************************************************************************
close all;

%Extract relevant parameters from raw traces
nFrames = size(rawTraces, 1);

%Enter user controlled visualization of individual particles 
t = 1;  %time
k = 2;  %Number of clusters

%Generate plots
[clusterFig, traceFig] = plotClusters(k, clusters{k-1}, rawTraces, time);
particleFig = plotFrame(t, clusters{k-1}, rawTraces);

%Colororder, can be edited if you want to make lines different colors
%colors = [lines(7); 0 0 0];
%colors = colors([3 1 2 4 5 6 7 8], :);

%Main control loop
while true
    %Get user input; do nothing unless input key is valid (see above)
    %Set dummy value for input key to enter loop that waits for input
    key = [];
    while isempty(key) || all(key ~= '.,mnxsh')
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
    
    %Increase/decrease number of clusters
    elseif key == 'm'
        if k < length(clusters)+1
            k = k+1;
            [clusterFig, traceFig] = ...
                plotClusters(k, clusters{k-1}, rawTraces, time, ...
                clusterFig, traceFig);
            %Clear histogram if open
            if exist('posHistFig', 'var')
                close(posHistFig);
                clear posHistFig
            end
        else
            beep;
            warning(['Max ', num2str(k), ' clusters'])
        end
    elseif key == 'n'
       if k > 2
            k = k-1;
            [clusterFig, traceFig] = ...
                plotClusters(k, clusters{k-1}, rawTraces, time, ...
                clusterFig, traceFig);
            %Clear histogram if open
            if exist('posHistFig', 'var')
                close(posHistFig);
                clear posHistFig
            end
        else
            beep;
            warning(['Min ', num2str(k), ' clusters'])
       end
    
    %Show position histogram   
    elseif key == 'h'
        if ~exist('posHistFig', 'var')
            posHistFig = ...
                figure('Units', 'Normalized', 'Position', [0 0 1, 0.45]);
            hold on;
            
            legendEntries = cell(1,k);
            for i = 1:k
                posHist = histcounts(nanmean(...
                    rawTraces(:,clusters{k-1}==i,3)), 0:0.01:1);
                plot(0.005:0.01:0.995, posHist,...
                    'LineWidth', 2);
                legendEntries{i} = ['Cluster ', num2str(i)];
            end
            
            title(['Position Histogram for k = ',num2str(k)])
            xlabel('AP Position')
            ylabel('Frequency')
            legend(legendEntries)
        else
            close(posHistFig)
            clear posHistFig
        end
        
    %Save movie
    elseif key == 's'
        %Initialize VideoWriter object
        filename = input('Enter Filename for Movie\n>>> ','s');
        label = input(['Enter title of movie\n',...
            '(Optional: to be displayed in each frame)\n>>> '], 's');
        
        %Create necessary folder
        if ~exist('Movies', 'dir')
            mkdir('Movies')
        end
        
        vidObj = VideoWriter(['Movies', filesep, filename]);
        vidObj.FrameRate = 8;
        open(vidObj);

        for tMovie = 1:nFrames
            %Generate figure
            particleFig = plotFrame(...
                tMovie, clusters{k-1}, rawTraces, particleFig, label);
            
            %Capture and write to movie
            frame = getframe(particleFig);
            writeVideo(vidObj, frame);
        end
        
        close(vidObj);
        
        figure(clusterFig);
        figure(traceFig);
        figure(particleFig);
    end
    
    particleFig = plotFrame(t, clusters{k-1}, rawTraces, particleFig);
end

close all

end

% Plotting Functions:

function particleFig = plotFrame(t,clusters, rawTraces, particleFig, label)
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
ax = gca;

for i = 1:nClusters
    whichParticles = (clusters == i);
    if ~isempty(whichParticles)
        plot(rawTraces(t, whichParticles, 3),...
            rawTraces(t, whichParticles, 5),'o');
    else
        ax.ColorOrderIndex = ax.ColorOrderIndex + 1;
    end
end

%Label with current time
timeStamp = text(0.21, 10, ['t = ', num2str(t-1), ' min']);
timeStamp.Color = fg;
timeStamp.FontWeight = 'bold';

%If an additional label is specified, include it in upper left corner
if nargin == 5
    labelText = text(0.21, 246, label);
    labelText.Color = fg;
    labelText.FontWeight = 'bold';
end

%Specify axis properties etc.
axis([0.2 0.85 0 256])
title(['Active nuclei colored by cluster; Frame: ', num2str(t)])
particleFig.Name = ['Frame ', num2str(t), ' of nc14'];
ax.Color = bg; ax.XColor = fg; ax.YColor = 'none';

particleFig = gcf;
end

function [clusterFig, traceFig] = ...
    plotClusters(k, clusters, rawTraces, time, clusterFig, traceFig)
%Make plots showing behavior of clusters over time
if nargin < 6
    clusterFig = figure('Units', 'Normalized', 'Position', [0 0 0.5 0.4]);
    traceFig = figure('Units', 'Normalized', 'Position', [0.5 0 0.5 0.4]);
end

%Plot each trace by cluster assignment
figure(clusterFig); cla;
[~, order] = sort(clusters);
imagesc(rawTraces(:,order,1)); colormap jet;
title(['Individual Traces: ' num2str(k), ' clusters'])
%Add lines separating clusters and labels
for i = 1:k
        hold on
        %Plot lines between clusters
        try
            if i ~= k
                plot(find(sort(clusters)==i, 1, 'last')*ones(1,2),...
                    [0 100], 'm', 'Linewidth', 2)
            end
            %Label each cluster
            text(mean(find(sort(clusters)==i)), 5, num2str(i), 'Color',...
                [1 1 1], 'FontWeight', 'bold', 'EdgeColor', [1 1 1]);
        catch
            % In case clusters aren't represented
        end
        
        ax = gca;
        ax.XTick = [];
        ax.XTickLabel = [];
        
        xlabel('Cluster')
        ylabel('Time')
end

%Plot mean traces for each cluster
figure(traceFig); cla
leg = cell(1,k);
fluoTrace = rawTraces(:,:,1);
fluoTrace(isnan(fluoTrace)) = 0;
for i = 1:k
    hold on
    plot(time, mean(fluoTrace(:,clusters == i),2),...
        'LineWidth', 2);
    xlabel('Frame Number')
    ylabel('Mean Fluorescence')
    leg{i} = ['Cluster ', num2str(i)];
end
legend(leg, 'Location', 'Best')
end