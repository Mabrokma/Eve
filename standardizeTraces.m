function [nc14, bin14, nc13, bin13] = ...
    standardizeTraces(genotype, bins)
%**************************************************************************
%
% Create standardized traces from the raw - query the linear interpolation
% to get flourescence etc. on one minute intervals (allow interval to be an
% input argument?).  Bin by 1% AP (also input argument?)
%
% Adapted from original extractTraces.m so that binning could be called w/o
% the entire extractTraces function
%
% This implementation is really stupid right now but it works
%
% Dependencies: none
% RW 8/2015
%
%**************************************************************************

if nargin ~= 2
    bins = genotype.CP.APbinID;
end

%Filter raw data for nc14--------------------------------------------------

queryTimes = 0:60;

rawTraces = genotype.rawTraces(genotype.CP.nc14:end,:,:);
filter14 = any(~isnan(rawTraces));
rawTraces = rawTraces(:,filter14(1,:,1),:);

dataTimes = genotype.CP.ElapsedTime ...
    - genotype.CP.ElapsedTime(genotype.CP.nc14);
dataTimes = dataTimes(dataTimes >= 0);

[nc14, bin14] = standardize(dataTimes, queryTimes, bins, rawTraces);

%Filter raw data for nc13--------------------------------------------------
if genotype.CP.nc13 ~= 0
    queryTimes = 0:25;
    
    rawTraces = ...
        genotype.rawTraces(genotype.CP.nc13:genotype.CP.nc14-1,:,:);
    filter13 = any(~isnan(rawTraces));
    rawTraces = rawTraces(:,filter13(1,:,1),:);
    
    dataTimes = genotype.CP.ElapsedTime ...
        - genotype.CP.ElapsedTime(genotype.CP.nc13);
    dataTimes = dataTimes(dataTimes >= 0 & ...
        dataTimes < dataTimes(genotype.CP.nc14));
    
    [nc13, bin13] = standardize(dataTimes, queryTimes, bins, rawTraces);
else
    nc13 = [];
    bin13 = [];
end

end

function [allTraces, binTraces] = ...
    standardize(dataTimes, queryTimes, bins, rawTraces)
%Remove NaNs in fluo
rawTraces(isnan(rawTraces(:,:,1))) = 0;
%Interpolate
allTraces = interp1(dataTimes, rawTraces, queryTimes);

%Add integrated fluorescence
allTraces(:,:,6) = cumsum(allTraces(:,:,1));
allTraces(:,:,7) = integrateWithDegradation(allTraces(:,:,1), 6);


%Bin traces to create a single 60 x 100 x 4 array
binTraces = zeros(length(queryTimes),length(bins)-1,3);

%**THIS METHOD KEEPS THE AP POSITION CONSTANT OVER TIME
% [~,~, whichBin] = histcounts(nanmean(allTraces(:,:,4)), bins);
% 
% for x = 1:length(bins)-1
%     %Total Fluorescence in bin
%     binTraces(:,x,1) = sum(allTraces(:,whichBin==x,1),2);
%     
%     %Mean nonzero fluorescence in bin
%     binTraces(:,x,2) = nanmean(allTraces(:,whichBin==x,1),2);
%     
%     %Total mRNA in bin at each time point
%     binTraces(:,x,3) = sum(allTraces(:,whichBin==x,7),2);
% end

%**THIS METHOD ALLOWS NUCLEI TO WIGGLE, IS SLOWER
for t = 1:length(queryTimes)
    [~,~, whichBin] = histcounts(allTraces(t,:,4), bins);
    
    for x = 1:length(bins)-1
        %Total Fluorescence in bin
        binTraces(t,x,1) = nansum(allTraces(t,whichBin==x,1));
        
        %Mean nonzero fluorescence in bin
        binTraces(t,x,2) = nanmean(allTraces(t,whichBin==x,1));
        
        %Total mRNA in bin at each time point
        binTraces(t,x,3) = nansum(allTraces(t,whichBin==x,7));
    end
end

end