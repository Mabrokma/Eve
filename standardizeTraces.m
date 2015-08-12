function [standardTraces, binTraces] = ...
    standardizeTraces(genotype, bins, queryTimes)
%**************************************************************************
%
% Create standardized traces from the raw - query the linear interpolation
% to get flourescence etc. on one minute intervals (allow interval to be an
% input argument?).  Bin by 1% AP (also input argument?)
%
% Adapted from original extractTraces.m so that binning could be called w/o
% the entire extractTraces function
%
% Dependencies: none (input must be a genotype structure from loadGenotype)
% RW 8/2015
%
%**************************************************************************

%Parse args
if nargin == 1
    bins = genotype.CP.APbinID;
    queryTimes = 0:60;
elseif nargin == 2
    queryTimes = 0:60;
end

%Filter raw data for nc14
rawTraces = genotype.rawTraces(genotype.CP.nc14:end,:,:);
filter14 = any(~isnan(rawTraces));
rawTraces = rawTraces(:,filter14(1,:,1),:);
rawTraces(isnan(rawTraces(:,:,1))) = 0;

dataTimes = genotype.CP.ElapsedTime ...
    - genotype.CP.ElapsedTime(genotype.CP.nc14);
dataTimes = dataTimes(dataTimes >= 0);

standardTraces = interp1(dataTimes, rawTraces, queryTimes);

%Add integrated fluorescence 
standardTraces(:,:,5) = cumsum(standardTraces(:,:,1));

%Bin traces to create a single 60 x 100 x 4 array
binTraces = NaN(length(queryTimes),length(bins)-1,3);
[~,~, bin] = histcounts(nanmean(standardTraces(:,:,4)),bins);

for x = 1:100
    %Total Fluorescence in bin
    binTraces(:,x,1) = sum(standardTraces(:,bin==x,1),2);
    
    %Mean nonzero fluorescence in bin
    binTraces(:,x,2) = nanmean(standardTraces(:,bin==x,1),2);
    
    %Total mRNA in bin at each time point
    binTraces(:,x,3) = sum(standardTraces(:,bin==x,5),2);
end