%Save heatmaps as individual figures
function [] = saveHeatmaps(genotype)
%Plot each type of heatmap on a different plot

nEmbryos = length(genotype);

APRange = [20 80];

titles = {'Sum Fluorescence in Bin', 'Mean Fluorescence of On Nuclei',...
    'Total mRNA Produced by Bin'};
saveName = {'SumFluo', 'MeanFluo', 'TotalmRNA'};

CLIMS = {[0 4E4] [0 4E3] [0 10E5]};

%Also make expression profiles- 
timePts = [10 20 30 40 50 59];
leg = cell(1,length(timePts));
for i = 1:length(timePts)
    leg{i} = ['t = ', num2str(timePts(i))];
end

for i = 1:nEmbryos
    for j = 1:2
        fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]);
        imagesc(APRange, [0 59], genotype(i).binTraces(...
            1:60,APRange(1):APRange(2),j), CLIMS{j});
        
        colormap jet;
        xlabel('AP Position (%EL)')
        ylabel('Time (min)')
        
        title(titles{j})
        
        saveas(fig, ['StripeFigures/HeatMaps/',...
            genotype(i).desc, num2str(i), '_', saveName{j}, '.png'])
        close(fig)
    end
    
    %Expression profile at each time point + composite
    for j = 1:6
        fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]);
        plot(APRange(1):APRange(2), ...
            genotype(i).binTraces(timePts(j),APRange(1):APRange(2),3),...
            'linewidth', 2);
        
        axis([APRange CLIMS{3}])
        
        xlabel('AP Position (%EL)')
        ylabel('Cumulative mRNA Produced (integrated Fluorescence)')
        title([titles{3}, '; t = ', num2str(timePts(j)), ' min'])
        
        saveas(fig, ['StripeFigures/HeatMaps/', genotype(i).desc,...
            't',num2str(timePts(j)), num2str(i), '_', saveName{3}, '.png'])
        close(fig)
    end
    
    fig = figure('Units', 'Normalized', 'Position', [0 0 1 1]);
    plot(APRange(1):APRange(2), ...
        genotype(i).binTraces(timePts,APRange(1):APRange(2),3)',...
            'linewidth', 2);
        
    axis([APRange CLIMS{3}])
    legend(leg);

    xlabel('AP Position (%EL)')
    ylabel('Cumulative mRNA Produced (integrated Fluorescence)')
    title(titles{3})

    saveas(fig,['StripeFigures/HeatMaps/', genotype(i).desc,...
        num2str(i), '_', saveName{3}, '.png'])
    close(fig)
end