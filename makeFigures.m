%make figures after loading all the necessary stuff

%%make individual heatmaps
close all
for i = 1:3
    figure; colormap jet
    imagesc([25, 80], [0, length(miniGFP(i).meanFluo)],...
        miniGFP(i).meanFluo(:,25:80,1)); colorbar;
    %Figure Properties
    title(['GFP - Embryo ', num2str(i), ': Mean Fluorescence'])
    xlabel('AP Position')
    ylabel('Frame Number')
    saveas(gca, ['StripeFigures/meanFluo_GFP', num2str(i), '.png']);
    close
    
    figure; colormap jet
    imagesc([25, 80], [0, length(miniGFP(i).meanFluo)],...
        miniGFP(i).meanFluo(:,25:80,2)); colorbar;
    %Figure Properties
    title(['GFP - Embryo ', num2str(i), ': Std Fluorescence'])
    xlabel('AP Position')
    ylabel('Frame Number')
    saveas(gca, ['StripeFigures/stdFluo_GFP', num2str(i), '.png']);
    close
    
    figure; colormap jet
    imagesc([25, 80], [0, length(miniGFP(i).meanFluo)],...
        miniGFP(i).activeNuclei(:,25:80,1),[0 1]); colorbar;
    %Figure Properties
    title(['GFP - Embryo ', num2str(i), ': Fraction of Active Nuclei ',...
        '(excluding disapproved)'])
    xlabel('AP Position')
    ylabel('Frame Number')
    saveas(gca, ['StripeFigures/activeNuc_GFP', num2str(i), '.png']);
    close
    
    figure; colormap jet
    imagesc([25, 80], [0, length(miniGFP(i).meanFluo)],...
        miniGFP(i).activeNuclei(:,25:80,2),[0 1]); colorbar;
    %Figure Properties
    title(['GFP - Embryo ', num2str(i), ': Fraction of Active Nuclei ',...
        '(including disapproved)'])
    xlabel('AP Position')
    ylabel('Frame Number')
    saveas(gca, ['StripeFigures/activeNucWD_GFP', num2str(i), '.png']);
    close

    figure; colormap jet
    imagesc([25, 80], [0, length(mininoGFP(i).meanFluo)],...
        mininoGFP(i).meanFluo(:,25:80,1)); colorbar;
    %Figure Properties
    title(['noGFP - Embryo ', num2str(i), ': Mean Fluorescence'])
    xlabel('AP Position')
    ylabel('Frame Number')
    saveas(gca, ['StripeFigures/meanFluo_noGFP', num2str(i), '.png']);
    close
    
    figure; colormap jet
    imagesc([25, 80], [0, length(mininoGFP(i).meanFluo)],...
        mininoGFP(i).meanFluo(:,25:80,2)); colorbar;
    %Figure Properties
    title(['noGFP - Embryo ', num2str(i), ': Std Fluorescence'])
    xlabel('AP Position')
    ylabel('Frame Number')
    saveas(gca, ['StripeFigures/stdFluo_noGFP', num2str(i), '.png']);
    close
    
    figure; colormap jet
    imagesc([25, 80], [0, length(mininoGFP(i).meanFluo)],...
        mininoGFP(i).activeNuclei(:,25:80,1), [0 1]); colorbar;
    %Figure Properties
    title(['noGFP - Embryo ', num2str(i), ': Fraction of Active Nuclei ',...
        '(excluding disapproved)'])
    xlabel('AP Position')
    ylabel('Frame Number')
    saveas(gca, ['StripeFigures/activeNuc_noGFP', num2str(i), '.png']);
    close
    
    figure; colormap jet
    imagesc([25, 80], [0, length(mininoGFP(i).meanFluo)],...
        mininoGFP(i).activeNuclei(:,25:80,2), [0 1]); colorbar;
    %Figure Properties
    title(['noGFP - Embryo ', num2str(i), ': Fraction of Active Nuclei ',...
        '(including disapproved)'])
    xlabel('AP Position')
    ylabel('Frame Number')
    saveas(gca, ['StripeFigures/activeNucWD_noGFP', num2str(i), '.png']);
    close
end

%average together mini cp structures and make heatmaps
%GFP-----------------------------------------------------------------------
GFPavg = averageMiniCP(miniGFP);
noGFPavg = averageMiniCP(mininoGFP);

figure; colormap jet
imagesc([25, 80], [0, length(GFPavg.meanFluo)],...
    GFPavg.meanFluo(:,25:80,1)); colorbar;
%Figure Properties
title('GFP - Consensus Trace: Mean Fluorescence')
xlabel('AP Position')
ylabel('Frame Number')
saveas(gca, 'StripeFigures/meanFluo_GFP_Consensus.png');
close

figure; colormap jet
imagesc([25, 80], [0, length(GFPavg.activeNuclei)],...
    GFPavg.activeNuclei(:,25:80,1), [0 1]); colorbar;
%Figure Properties
title(['GFP - Consensus Trace: Fraction of Active Nuclei ',...
    '(excluding disapproved)'])
xlabel('AP Position')
ylabel('Frame Number')
saveas(gca, 'StripeFigures/activeNuc_GFP_Consensus.png');
close

%noGFP---------------------------------------------------------------------
figure; colormap jet
imagesc([25, 80], [0, length(noGFPavg.meanFluo)],...
    noGFPavg.meanFluo(:,25:80,1)); colorbar;
%Figure Properties
title('noGFP - Consensus Trace: Mean Fluorescence')
xlabel('AP Position')
ylabel('Frame Number')
saveas(gca, 'StripeFigures/meanFluo_noGFP_Consensus.png');
close

figure; colormap jet
imagesc([25, 80], [0, length(noGFPavg.activeNuclei)],...
    noGFPavg.activeNuclei(:,25:80,1), [0 1]); colorbar;
%Figure Properties
title(['noGFP - Consensus Trace: Fraction of Active Nuclei ',...
    '(excluding disapproved)'])
xlabel('AP Position')
ylabel('Frame Number')
saveas(gca, 'StripeFigures/activeNuc_noGFP_Consensus.png');
close

%subtract gfp and nogfp average heatmaps to find differences---------------
traceLength = min(length(GFPavg.meanFluo), length(noGFPavg.meanFluo));

meanFluoDiff = GFPavg.meanFluo(1:traceLength,:,:) ...
    - noGFPavg.meanFluo(1:traceLength,:,:);
activeNucleiDiff = GFPavg.activeNuclei(1:traceLength,:,:) ...
    - noGFPavg.activeNuclei(1:traceLength,:,:);

figure; colormap jet(5)
imagesc([25, 80], [0, length(meanFluoDiff)],...
    meanFluoDiff(:,25:80,1)); colorbar;
%Figure Properties
title('Differences: Mean Fluorescence')
xlabel('AP Position')
ylabel('Frame Number')
saveas(gca, 'StripeFigures/meanFluo_DIFF.png');
close

figure; colormap jet(5)
imagesc([25, 80], [0, length(activeNucleiDiff)],...
    activeNucleiDiff(:,25:80,1), [0 1]); colorbar;
%Figure Properties
title(['Differences: Fraction of Active Nuclei ',...
    '(excluding disapproved)'])
xlabel('AP Position')
ylabel('Frame Number')
saveas(gca, 'StripeFigures/activeNuc_DIFF.png');
close

%cluster traces