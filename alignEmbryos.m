function [alignedGenotype, avgAligned, allAligned] = ...
    alignEmbryos(genotype, objfun)
%**************************************************************************
% Take an entire genotype of processed data and slightly shift or scale to
% minimize the covariance between binned expression patterns (over time).
% The logic is that there is sufficient error in the viewing window vs.
% full embryo and in the variance of the pattern itself, so if we want to
% consider stripe-specific or stripe-border-specific qualities of traces,
% we want respective stripe-n's to agree as much as possible in expression.
%
% We limit alignment to minimal shifts (+/- 2%AP) and scales (+/-5%)
%
% The optimization is performed with simmulated annealing on a
% user-specified objective function, optionally passed as a second
% argument. The defualt is mean squared difference between images
% normalized to their respective maximum, included as a subfunction,
% however any function can be passed that has the syntax 
%   f(shift, scale, template, genotype)
% where template is the image to which the other embryos are aligned
%
% /Eve/Metrics/ contains a few other image comparison metrics, including
% ssim, corr2, and a binary agreement metric
% RW 8/2015
%**************************************************************************
totalTime = tic;

if nargin < 2
    objfun = @calculateMSD;
end

%Limits on shifting and scaling
maxShift = 0.02;
maxScale = 0.1;


%Explore parameter space using simulated annealing

%Use standardTraces (for nc14) to align.  Do correlation on 1% bins,
%repeatedly shift/scale until all combinations are exhausted, keep track of
%best alignment
template = genotype(1).binTraces(:,:,1);
alignedGenotype = genotype;

usedShifts = zeros(1,length(genotype)-1);
usedScales = zeros(1,length(genotype)-1);
function_vals = zeros(1,length(genotype)-1);

fprintf('\nAlgining embryos using "%s" as objective function...\n',...
    func2str(objfun));

parfor i = 2:length(genotype)
    thisTime = tic;
    %Create anonymous function with only shift and scale as parameters
    %Must be treated as a single vector for optimization w/simulannealbnd()
    objfun_anon = @(x) objfun(x(1), x(2), template, genotype(i));
    
    %Call matlab's simmulated annealing function, keep track of optimal
    %shifts and scales (params) and the value of the objective function, as
    %well as the output structure from simulannealbnd()
    [params, function_vals(i-1), ~, output] = simulannealbnd(...
        objfun_anon, [0 0], [-maxShift 1-maxScale], [maxShift 1+maxScale]);
    
    %Get optimal shifts and scales
    bestShift = params(1);
    bestScale = params(2);
    
    %Display
    fprintf('(Embryo %i: %i iterations in %4.1f seconds)\n',...
        i, output.iterations, toc(thisTime))

    alignedGenotype(i).rawTraces(:,:,3:4) = ...
        (2-bestScale)*genotype(i).rawTraces(:,:,3:4) - bestShift;
    
    [alignedGenotype(i).standardTraces, alignedGenotype(i).binTraces] = ...
        standardizeTraces(alignedGenotype(i));

    usedShifts(i-1) = bestShift;
    usedScales(i-1) = bestScale;
end

%Print results
fprintf('\nAlignment complete!\n\n')
fprintf('  %i Embryos aligned in %4.1f seconds\n',...
    length(genotype)-1, toc(totalTime));
fprintf('\n          Shift   | Scale  | Score\n')
for i = 1:length(genotype)-1
    fprintf('Embryo %i: % 4.4f | %5.4f | %4.3f\n',...
        i+1, usedShifts(i), usedScales(i), function_vals(i))
end
fprintf('\n')

%Generate figure to show how the alignment went
[avgAligned, allAligned] = averageTraces(alignedGenotype);
avgUnaligned = averageTraces(genotype);

scale = max([max(max(avgUnaligned(:,:,1)))...
    max(max(avgAligned(:,:,1)))]); %lol

figure('Units', 'Normalized', 'Position', [0 0 1 1]);
colormap jet;
subplot(1,2,1)
imagesc([20 85], [0 60], avgUnaligned(:,20:85,1), [0 scale]);
title('Unaligned average heatmap')

subplot(1,2,2)
imagesc([20 85], [0 60], avgAligned(:,20:85,1), [0 scale]);
title('Aligned average heatmap')

end

function msd = calculateMSD(shift, scale, template, genotype)
%Calculate image with given parameters,
[~, test] = standardizeTraces(genotype, scale*(0:0.01:1) + shift);
test = test(:,:,1);

%Compare to unshifted image
sd = (template/max(template(:)) - test/max(test(:))).^2;
msd = nanmean(sd(:));
end