function [] = plotHeatmaps(genotype)
%Plot each type of heatmap on a different plot

nEmbryos = length(genotype);
nRows = ceil(nEmbryos / 3);

APRange = [20 80];

titles = {'Sum Fluorescence in Bin', 'Mean Fluorescence of On Nuclei',...
    'Total mRNA Produced by Bin'};
CLIMS = {[0 4E4] [0 4E3] [0 10E5]};

for j = 1:3
    figure('Units', 'Normalized', 'Position', [0 0 1 1]);
    for i = 1:nEmbryos
        subplot(nRows,3,i)
        imagesc(APRange, [0 59], genotype(i).binTraces(...
            1:60,APRange(1):APRange(2),j), CLIMS{j});
        colormap jet;
        xlabel('AP Position')
        ylabel('Time')
        
        if i == 2
            title(titles{j})
        end
    end
    
end