function dif = embryo_logical(shift, scale, template, genotype)
%Calculate image with given parameters,
[~, test] = standardizeTraces(genotype, scale*(0:0.01:1) + shift);
test = test(:,:,1)~=0;
template = template ~= 0;

%Compare to unshifted image
sd = (template - test).^2;
dif = sum(sd(:));
end
