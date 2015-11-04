function [matches] = motifScan(rawData, motif, R, distance_fun)
%**************************************************************************
%
% Look for occurrences of a user-specified motif in a time-series.  A motif
% match is defined by a euclidean distance (or custom distance function) of
% less than R between a test subsequence and the supplied motif.
%
%**************************************************************************
tic;
timeSeries = rawData(:,:,1);
pos = nanmean(rawData(:,:,3));

%Check validity of input motif
l = length(motif);
traceLength = size(timeSeries,1);
totalLength = numel(timeSeries);

if any(isnan(motif)) || l == 0 || l > traceLength
    error('Invalid motif')
end

taxicab = @(x,y)sum(abs(minus(x,y)));
if nargin < 3
    %Default metric is taxicab
    distance_fun = taxicab;
    R = 2E2*l;
elseif nargin < 4
    distance_fun = taxicab;
end

emptyThreshold = min(0.6,sum(motif==0)/l);

%Initialize counter at max capacity
matches = NaN(floor(totalLength/l - traceLength*l),3);
hits = 0;

for i = 1:totalLength - l + 1
    
    %Subsequence i(l)
    iIdx = i:i+l-1;            %indices of subsequence i
    c_i = timeSeries(iIdx);    %values of subsequence i
    
    %Skip bad subsequences
    skipFlag = (any(mod(iIdx,traceLength) == 0) ...
        && any(mod(iIdx,traceLength) == 1)) ...  cross traces
        || sum(c_i == 0)/l > emptyThreshold;%    empty
    
    if skipFlag
        continue
    end
    
    d = distance_fun(c_i,motif);
    if d <= R
        %Count occurence of motif
        hits = hits+1;
        matches(hits,1) = i;
        matches(hits,2) = d;
        if isequal(distance_fun,taxicab)
            matches(hits,3) = d;
        else
            matches(hits,3) = taxicab(c_i,motif);
        end
    end    
end

%Filter out mutliple overlapping matches
matches(isnan(matches(:,1)),:) = [];

elapsedTime = toc;
fprintf('Scanning completed in %5.2f seconds\n', elapsedTime)
%Plot results--------------------------------------------------------------

figure('Units', 'Normalized', 'Position', [0 0 1 1]);
subplot(2,2,1:2)
for i = 1:length(matches)
    hold on
    try
        plot(0:l+1,timeSeries(matches(i,1)-1:matches(i,1)+l))
    catch
        plot(1:l,timeSeries(matches(i,1):matches(i,1)+l-1))
    end
end
plot(motif, 'k', 'linewidth', 3);
title(['Motif - ', num2str(length(matches)), ' Matches'])
xlabel('Relative time')
ylabel('Fluorescence')
plot([1,1],[0,1.2*max(motif)],'k--');
plot([l,l],[0,1.2*max(motif)],'k--');
axis([0, l+1, 0, 1.2*max(motif)+1])

%Heatmap of motif occurance
subplot(2,2,3)
bins = 40;
motifMap = NaN(traceLength,bins);
mods = mod(matches(:,1), traceLength);
pos = pos(ceil(matches(:,1)/traceLength));
for t = 1:traceLength
    if t ~= traceLength
        whichPos = pos(mods==t);
    else
        whichPos = pos(mods==0);
    end
    
    %fprintf('%i Occurences at t = %i\n', length(whichPos), t-1)

    for x = 1:bins
        motifMap(t,x) = sum(whichPos > (x-1)/bins & whichPos <= x/bins);
    end
end
imagesc([0 100], [0 60], motifMap);
colormap jet; colorbar;
title('Positional and Spatial occurance of Motif')
xlabel('AP Position (%EL)')
ylabel('Time')

%Histogram of distances
subplot(2,2,4)
histogram(matches(:,3));
title('Histogram of distances from motif')
xlabel('Cityblock distance')
ylabel('Frequency')