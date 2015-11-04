function [whichCluster, centroids, fluoTrace, time] = ...
    clusterTraces(fluoTrace, time, varargin)
%**************************************************************************
% (Under development)
%
% Takes a (single) genotype struct and clusters individual fluorescence
% traces using kmeans with k = 2:8, makes a bunch of figures.  nucOn14[] is
% a logical vector specifying on/off times, can be used in clustering?
% 
% Legacy support for genotype structure may not be necessary / functional
%
% Dependencies: none
% RW 8/2015
%**************************************************************************

%Allow for entering of a genotype structure instead of a trace/time combo
if nargin == 1
    genotype = fluoTrace;
    %Extract nuclear cycle 14 from raw traces
    fluoTrace= genotype.rawTraces(genotype.CP.nc14:genotype.CP.nc14+100,...
        [genotype.CP.CompiledParticles.nc] == 14,:);
    time= genotype.CP.ElapsedTime(genotype.CP.nc14:genotype.CP.nc14+100)...
        - genotype.CP.ElapsedTime(genotype.CP.nc14);
end


%Control generation of figures
plotFlag = 1;
saveFlag = 0;
metric = 'sqeuclidean';
validMetrics = ...
    {'sqeuclidean', 'cityblock', 'correlation', 'cosine', 'hamming'};
if nargin > 2
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'noplot')
            plotFlag = 0;
        elseif strcmpi(varargin{i}, 'save')
            saveFlag = 1;
            saveFolder = input('Name save folder\n>?','s');
            if ~exist(['Figures/Clustering/', saveFolder], 'dir')
                mkdir(['Figures/Clustering/', saveFolder])
            else
                warning(['Folder already exists,',...
                    'this may overwrite existing figures.'])
                fprintf('Press any key to continue...\n')
                pause;
            end
        elseif any(cellfun(@(x)strcmpi(x,varargin{i}), validMetrics))
            metric = varargin{i};
        end
    end
end

if plotFlag
    %Histogram nonzero fluorescences
    binsize = 100;
    fluoVector = fluoTrace(:,:,1);
    [fluoHist, edges] = histcounts(fluoVector(:), binsize);
    centers = edges + binsize/2; centers(end) = [];
    figure;
    plot(centers,fluoHist)
    axis([0 12000 0 750])
    title('Histogram of (nonzero) fluorescence intensities')
    xlabel('Fluo')
    ylabel('Frequency')
    if saveFlag
        export_fig(sprintf('Figures/Clustering/%s/hist.png', saveFolder))
    end
end


if strcmpi(metric, 'hamming')
    fluoTrace = ~isnan(fluoTrace);
else
    fluoTrace(isnan(fluoTrace(:,:,1))) = 0;
end

colors = [lines(7); 0 0 0];
maxClusters = 8;
whichCluster = cell(1, maxClusters-1);
centroids = cell(1, maxClusters-1);

h = waitbar(0);
for k = 2:maxClusters
    waitbar((k-1)/(maxClusters-1), h, ['Clustering with k = ', num2str(k)])
    %Cluster on the individual traces   
    [whichCluster{k-1}, centroids{k-1}] = kmeans(fluoTrace(:,:,1)',...
        k, 'Replicates', 10, 'distance', metric);
    
    [~, order] = sort(whichCluster{k-1});
    
    if plotFlag
        figure('Units', 'Normalized', 'Position', [0 0 1 1]);
        subplot(2,2,1)
        imagesc(fluoTrace(:,order,1)); colormap jet
        title(['Individual Traces: ' num2str(k), ' clusters'])
        for i = 1:k
            hold on
            %Label each cluster
            text(mean(find(sort(whichCluster{k-1})==i)), 5, num2str(i),...
                'Color', [1 1 1], 'FontWeight', 'bold',...
                'EdgeColor', [1 1 1]);
            
            %Plot lines between clusters
            if i ~= k
                plot(find(sort(whichCluster{k-1})==i, 1, 'last')...
                    *ones(1,2), [0 100], 'm', 'Linewidth', 2)
            end
            ax = gca;
            ax.XTick = [];
            ax.XTickLabel = [];
            
            xlabel('Clusters')
            ylabel('Time')
        end
        
        %Plot mean traces for each cluster & bin clusters by AP Position
        leg = cell(1,k);
        for i = 1:k
            subplot(2,2,2)
            hold on
            plot(time, nanmean(fluoTrace(:,whichCluster{k-1} == i),2),...
                'Color',colors(i,:), 'LineWidth', 2);
            xlabel('Time')
            ylabel('Mean Fluorescence')
            leg{i} = ['Cluster ', num2str(i)];
            
            subplot(2,2,[3 4])
            clusterHist = histcounts(nanmean(...
                fluoTrace(:,whichCluster{k-1}==i,3)), 0:0.01:1);
            hold on
            plot(0.005:0.01:0.995, clusterHist,...
                'Color', colors(i,:), 'LineWidth', 2)
            xlabel('AP Position')
            ylabel('Frequency')
        end
        subplot(2,2,2)
        title('Mean traces by cluster')
        legend(leg, 'Location', 'Best')
        subplot(2,2,[3 4])
        title('Position histogram by Cluster')
        legend(leg, 'Location', 'Best')
        if saveFlag
            export_fig(sprintf('Figures/Clustering/%s/k%i.png',...
                saveFolder, k))
        end
    end
end
delete(h)


