function [handle] = compareExpression(varargin)
%Compare wt, kni/+, and kni/kni over time

timepts = [5 15 25 35 45 55];
leg = cell(1,length(timepts));

for i = 1:length(timepts)
    leg{i} = ['t = ', num2str(timepts(i))];
end

AP = [15 85];
RANGE1 = [0 14E5];
RANGE2 = [0 10E5];

handle = figure('Units', 'Normalized', 'Position', [0 0 1 1]);
colormap jet;

subplot(2,3,4)
plot(dorsal(2).binTraces(timepts,:,3)')
title('wt, dorsal')
xlabel('AP Position (%El)')
ylabel('Cumulative mRNA (AU)')
axis([AP RANGE1])
legend(leg, 'Location', 'Northeast')

subplot(2,3,1)
imagesc(AP, 0:59, dorsal(2).binTraces(:,AP(1):AP(2),3), RANGE2);
%colorbar
title('wt, dorsal')
xlabel('AP Position (%El)')
ylabel('Time (min)')

subplot(2,3,5)
plot(avgKniWt(timepts,:,3)')
title('kni/+, dorsal')
xlabel('AP Position (%El)')
ylabel('Cumulative mRNA (AU)')
axis([AP RANGE1])
legend(leg, 'Location', 'Northeast')


subplot(2,3,2)
imagesc(AP, 0:59, avgKniWt(:,AP(1):AP(2),3), RANGE2);
%colorbar;
title('kni/+, dorsal')
xlabel('AP Position (%El)')
ylabel('Time (min)')

subplot(2,3,6)
plot(avgKniKni(timepts,:,3)')
title('kni/kni, dorsal')
xlabel('AP Position (%El)')
ylabel('Cumulative mRNA (AU)')
axis([AP RANGE1])
legend(leg, 'Location', 'Northeast')

subplot(2,3,3)
imagesc(AP, 0:59, avgKniKni(:,AP(1):AP(2),3), RANGE2);
%colorbar('SouthOutside');
title('kni/kni, dorsal')
xlabel('AP Position (%El)')
ylabel('Time (min)')

end

function [] = plotHeatmap(data, theTitle)
    imagesc(AP, 0:59, data, RANGE2);
    %colorbar('SouthOutside')
    title(theTitle)
    xlabel('AP Position (%El)')
    ylabel('Time (min)')
end

function [] = plotContour(data, theTitle)
    plot(data)
    title(theTitle)
    xlabel('AP Position (%El)')
    ylabel('Cumulative mRNA (AU)')
    axis([AP RANGE1])
    legend(leg, 'Location', 'Northeast')
end
