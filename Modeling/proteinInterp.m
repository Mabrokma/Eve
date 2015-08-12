%Interpolate a set of protein concentrations given at discrete times to
%produce a value for each time in ElapsedTime
%inProteins should be a t x 100 x n matrix with t time steps, in 100 ap
%bins and n different proteins

function [linearInterp] = ...
    proteinInterp(ElapsedTime, tfConc, proteinTimes)

%Do linear
linearInterp = ...
    interp1(proteinTimes, tfConc, ElapsedTime, 'linear');

figure(1)
subplot(1,2,1)
imagesc(tfConc)
title('Input Concentrations')
subplot(1,2,2)
imagesc(linearInterp(1:find(ElapsedTime > proteinTimes(end),1),:))
title('Linear Interpolation')
colormap jet; colorbar