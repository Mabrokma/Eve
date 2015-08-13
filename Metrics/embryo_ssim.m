function ssimval = embryo_ssim(shift, scale, template, genotype)
%Calculate image with given parameters,
[~, test] = standardizeTraces(genotype, scale*(0:0.01:1) + shift);
test = test(:,:,1);

%Compare to unshifted image
ssimval = ssim(test, template);
end