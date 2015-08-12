function [whichCluster, nc14, time14] = clusterTraces(genotype, varargin)
%**************************************************************************
% (Under development)
%
% Takes a (single) genotype struct and clusters individual fluorescence
% traces using kmeans with k = 2:8, makes a bunch of figures.  Metric is
% specified in the kmeans command (line 44ish) and can be changed,
% currently sqeuclidean.  nucOn14[] is a binary vector specifying on/off
% times, can be used in clustering?
% 
% Currently a rough version that isn't super readable, I intend to clean up
% all the multi-level struct references
%
% Dependencies: none
% RW 8/2015
%**************************************************************************

%Control generation of figures
figureFlag = 1;
metric = 'sqeuclidean';
if nargin > 1
    for i = 1:length(varargin)
        if strcmpi(varargin{i}, 'nofigures')
            figureFlag = 0;
        end
    end
end

if figureFlag
    %Histogram nonzero fluorescences
    fluoVector = genotype.rawTraces(:,:,1); fluoVector = fluoVector(:);
    binsize = 100;
    [fluoHist, edges] = histcounts(fluoVector, binsize);
    centers = edges + binsize/2; centers(end) = [];
    figure;
    plot(centers,fluoHist)
    axis([0 12000 0 750])
    title(['Histogram of (nonzero) fluorescence intensities: ',...
        genotype.Name])
end

%Extract nuclear cylce 14 from raw traces
nc14 = genotype.rawTraces(genotype.CP.nc14:genotype.CP.nc14+100,...
    [genotype.CP.CompiledParticles.nc] == 14,:);
fluonc14 = nc14(:,:,1)';
time14 = genotype.CP.ElapsedTime(genotype.CP.nc14:genotype.CP.nc14+100) ...
    -genotype.CP.ElapsedTime(genotype.CP.nc14);

%nucOn14 = ~isnan(fluonc14); %Logical array of on/off time points
fluonc14(isnan(fluonc14)) = 0;

colors = [lines(7); 0 0 0];
maxClusters = 8;
whichCluster = cell(1, maxClusters);

h = waitbar(0);
for k = 2:maxClusters
    waitbar((k-1)/(maxClusters-1), h, ['Clustering with k = ', num2str(k)])
    %Cluster on the individual traces   
    whichCluster{k-1} = kmeans(fluonc14, k, 'Replicates', 10, ...
        'distance', metric);
    [~, order] = sort(whichCluster{k-1});
    
    if figureFlag
        figure('Units', 'Normalized', 'Position', [0 0 1 1]);
        subplot(2,2,1)
        imagesc(fluonc14(order,:)'); colormap jet
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
            plot(time14, nanmean(fluonc14(whichCluster{k-1} == i, :)),...
                'Color',colors(i,:), 'LineWidth', 2);
            xlabel('Time')
            ylabel('Mean Fluorescence')
            leg{i} = ['Cluster ', num2str(i)];
            
            subplot(2,2,[3 4])
            clusterHist = histcounts(...
                nanmean(nc14(:,whichCluster{k-1}==i,3)), 0.2:0.01:0.85);
            hold on
            plot(0.205:0.01:0.845, clusterHist,...
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
    end
end
delete(h)

%** TODO 

%Plot all individual traces in colors corresponding to stripe assignment??

%Plot all individual traces in colors corresponding to location in stripe??