function [normDifExp, difExp] = compareGapGenes(gapExp, t, x)
%Compare the expression profiles of the gap genes over time and space to a
%given point

if nargin < 2
    t = 35;
end
 
if nargin < 3
    x = [33 41 46 52 59 69 78];
end

nPts = length(x);
gapSize = size(gapExp);

difExp = NaN(gapSize(1), gapSize(2), nPts);
normDifExp = NaN(gapSize(1), gapSize(2), nPts);

for i = 1:nPts
    ptExp = repmat(gapExp(t,x(i),:), gapSize(1), gapSize(2));
    difExp(:,:,i) = sum(abs(gapExp - ptExp),3);
    normDifExp(:,:,i) = difExp(:,:,i) ./ sum(gapExp,3);
    
%     %Plot results
%     figure('Units', 'Normalized', 'Position', [0 0 1 1]);
%     imagesc(normDifExp(:,:,i), [0 2]);
%     colormap jet; colorbar;
%     hold on;
%     
%     plot(x(i), t, 'xw');
%     title(['Similarity of gap gene profiles to Stripe ', ...
%         num2str(i), ' (point "x")'])
    %saveas(gca,['GapGeneStripe', num2str(i), '.png'])
end

%Sum Reciprocal similarity
difOverlay = sum(1./(0.1+normDifExp),3);
figure('Units', 'Normalized', 'Position', [0 0 1 1])
imagesc(difOverlay)%, [0 10])
colormap jet; colorbar;
title(['Sum Recriprocal Similarity - Seed time: ', num2str(t), 'm'])
%saveas(gca, ['AllSumRecip_t', num2str(t), '.png'])
