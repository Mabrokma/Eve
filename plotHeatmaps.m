function [] = plotHeatmaps(genotype, varargin)
%Plot each type of heatmap on a different plot

saveFlag = 0;
%Parse additional args
if ~isempty(varargin)
    for i = 1:length(varargin)
        if strcmp(varargin{i}, 'save')
            saveFlag = 1;
            desc = genotype(1).desc;
            savePath = ['Figures',filesep,'Heatmaps',filesep,desc,filesep];
            if ~exist(savePath, 'dir')
                mkdir(savePath)
            end
        end
    end
end



%Get info from data
nEmbryos = length(genotype);
if nEmbryos <= 9
    nRows = ceil(nEmbryos / 3);
    individualPlots = 0;
else
    individualPlots = 1;
end

%Plotting features
titles = {'Sum Fluorescence in Bin', 'Mean Fluorescence of On Nuclei',...
    'mRNA Accumulated (t50 = 6min)'};
saveNames = {'sumFluo','meanFluo','mRNA'};
APRange = [20 80];
CLIMS = {[0 4E4], [0 5E3], [0 1E5]};

if individualPlots
    % Plot many embryos in different plots
    for i = 1:nEmbryos
        figure('Units', 'Normalized', 'Position', [0 0.25 1 0.5]);
        for j = 1:3
            subplot(1,3,j)
            imagesc(APRange, [0 59], genotype(i).bin14(...
                1:60,APRange(1):APRange(2),j));
            colormap jet; colorbar;
            title(titles{j})
            xlabel('AP Position')
            ylabel('Time')
        end
        if saveFlag
            export_fig([savePath, desc,'_',num2str(i),'.png']);
            close;
        end 
    end
else
    % If there are sufficiently few, plot them all together
    for j = 1:3
        figure('Units', 'Normalized', 'Position', [0 0 1 1]);
        for i = 1:nEmbryos
            subplot(nRows,3,i)
            imagesc(APRange, [0 59], genotype(i).bin14(...
                1:60,APRange(1):APRange(2),j), CLIMS{j});
            colormap jet;
            xlabel('AP Position')
            ylabel('Time')
            
            if i == 2
                title(titles{j})
            end
        end
        
        if saveFlag
            export_fig([savePath, desc,'_',saveNames{j},'.png']);
            close;
        end 
    end
end