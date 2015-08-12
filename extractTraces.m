function [rawTraces, standardTraces, binTraces, binAllTraces] = ...
    extractTraces(cp, nc)
%**************************************************************************
%Extract traces from a given compiled particles structure
%
%Everything (raw) is stored in a single array:
%   rawTraces(:,:,1) - Fluorescence
%   rawTraces(:,:,2) - standard deviation of fluorescence
%   rawTraces(:,:,3) - AP position at time point
%   rawTraces(:,:,4) - AP position of corresponding nucleus
%   rawTraces(:,:,5) - y position (in pixels) of particle in image
%                   (use it for displaying particles)
%
% rawTraces include every particle tracked and it's info
%
% standardTraces queries a linear interpolation of the rawTraces onto the
%    first 60 minutes of nc14
%
% binTraces bins the standardTraces by 1% AP
%
% binAllTraces bins rawTraces by 1% AP
% 
% Dependencies: particlesInFrame.m
% RW 7/2015
%**************************************************************************

%Pares nc arg: only affects binTraces/ standardTraces
if nargin == 1
    nc = 14;
end
ncName = ['nc', num2str(nc)];

nParticles = length(cp.CompiledParticles);
traceLength = length(cp.ElapsedTime);

rawTraces = NaN(traceLength, nParticles, 4);

%Get positions and fluorescences of each particle in each frame
for t = 1:length(cp.ElapsedTime)
    [idx,~,fluo,fluoError,pPos,nPos,yPos] = ...
        particlesInFrame(cp.CompiledParticles,t,...
        'Fluo','FluoError','APposParticle','APpos', 'yPos');
    
    if ~isempty(idx) 
        rawTraces(t, idx, 1) = fluo;
        rawTraces(t, idx, 2) = fluoError;
        rawTraces(t, idx, 3) = pPos;
        rawTraces(t, idx, 4) = nPos;
        rawTraces(t, idx, 5) = yPos;
    end
end

%Interpolate (single-frame) holes in the data
for i = 1:nParticles
    for t = 2:traceLength-1
        if isnan(rawTraces(t,i,1)) ...
                && ~isnan(rawTraces(t+1,i,1)) ...
                && ~isnan(rawTraces(t-1,i,1))
            rawTraces(t,i,:) = (rawTraces(t+1,i,:) + rawTraces(t-1,i,:))/2;
        end
    end
end

%Bin raw traces
binAllTraces = zeros(traceLength,100,2);
[~,~, bin] = histcounts(nanmean(rawTraces(:,:,4)),cp.APbinID);
for t = 1:length(cp.ElapsedTime)
    for x = 1:100
        %Total Fluorescence in bin
        binAllTraces(t,x,1) = nansum(rawTraces(t,bin==x,1));
        %Mean nonzero fluorescence in bin
        binAllTraces(t,x,2) = nanmean(rawTraces(t,bin==x,1));
    end
end

%Create the standardized trace from individual traces
%Janky way of looking at nc13 also (** TODO: Clean this up)
queryTimes{2} = 0:59;
queryTimes{1} = 0:15;

dataTimes = cp.ElapsedTime - cp.ElapsedTime(cp.(ncName));
if isfield(cp, ['nc', num2str(nc+1)])
    endFrame = cp.(['nc', num2str(nc+1)]);
else
    endFrame = length(cp.ElapsedTime);
end

rawTraces(isnan(rawTraces(:,:,1))) = 0;

standardTraces = interp1(dataTimes(cp.(ncName):endFrame), ...
    rawTraces(cp.(ncName):endFrame,:,:), queryTimes{nc-12});

%Add integrated fluorescence 
standardTraces(:,:,5) = cumsum(standardTraces(:,:,1));

rawTraces(rawTraces == 0) = NaN;

%Bin traces to create a single 60 x 100 x 4 array
binTraces = NaN(length(queryTimes{nc-12}),length(cp.APbinID)-1,3);
[~,~, bin] = histcounts(nanmean(standardTraces(:,:,4)),cp.APbinID);
%Bin nuclei for calculating fractional activity
if isfield(cp, 'EllipsesFilteredPos')
    [~,~, nucBin] = ...
        histcounts(nanmean(cp.EllipsesFilteredPos{t}), cp.APbinID);
else
    [~,~, nucBin] = ...
        histcounts(nanmean(cp.EllipsePos{t}), cp.APbinID);
end

for t = queryTimes{nc-12} + 1;
    for x = 1:100
        %Total Fluorescence in bin
        binTraces(t,x,1) = sum(standardTraces(t,bin==x,1));
        
        %Mean nonzero fluorescence in bin
        binTraces(t,x,2) = nanmean(standardTraces(t,bin==x,1));
        
        %Total mRNA in bin at point t
        binTraces(t,x,3) = sum(standardTraces(t,bin==x,5));
        
        %Fraction of active Nuclei
        nP = nansum(standardTraces(t,bin==x,5));
        nN = sum(nucBin == x);
        binTraces(t,x,4) = nP/nN;
    end
end
