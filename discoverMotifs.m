function [motifs] = discoverMotifs(timeSeries, l, R, distance_fun)
%**************************************************************************
%
% As an alternative to clustering, find maximally agreeing subsequences in
% a time-series, excluding trivial matches (i.e. shifts of one subsequence
% to each other).  Based loosely on the brute-force algorithm presented in
% Lin et al., 2002.
%
% Search for motifs of increasing length, checking every subsequence
% against every other non-overlapping subsequence and declaring the densest
% cluster to be the 1-motif, and jth dense the j-motif of the time-series.
%
% Distance function for comparing subsequences can be specified as an
% argument
%
% TODO - don't save all motifs, only top 5 (or something), dynamically
% check whether a given motif has enough instances to be in the top 5
%**************************************************************************

if nargin < 3
    %Default metric is euclidean distance
    distance_fun = @(x,y)sum(power(minus(x,y),2));
    R = 1E4;
elseif nargin < 4
    distance_fun = @(x,y)sum(power(minus(x,y),2));
end

displayFlag = 1;

%Limits on sequence lengths
traceLength = size(timeSeries,1);
totalLength = numel(timeSeries);

emptyThreshold = 0.6;

motifs = cell(1,totalLength-traceLength*l);

h = waitbar(0);
for i = 1:totalLength - l + 1
    %Initialize count & preallocate index matrix
    matches = 0;
    motifs{i} = NaN(1,floor(totalLength/l - traceLength*l));
    
    %Subsequence i(l)
    iIdx = i:i+l-1;            %indices of subsequence i
    c_i = timeSeries(iIdx);    %values of subsequence i
    
    %Skip bad subsequences
    skipFlag = (any(mod(iIdx,traceLength) == 0) ...
        && any(mod(iIdx,traceLength) == 1)) ...  cross traces
        || sum(c_i == 0)/l > emptyThreshold;%   empty
    
    if skipFlag
        motifs{i} = [];
        continue
    end
    
    for j = 1:totalLength - l + 1;
        waitbar(i*j / (totalLength-l + 1)^2, h, 'Finding motifs...')
        %Subsequence j(l)
        jIdx = j:j+l-1;
        c_j = timeSeries(jIdx);
        
        %Skip bad subsequences again
        skipFlag = (any(j == iIdx) || any(i == jIdx)) ...   overlapping
            || (any(mod(jIdx,traceLength) == 0) ...         
            && any(mod(jIdx,traceLength) == 1)) ...         cross traces
            || sum(c_j == 0)/l > emptyThreshold;%           empty
        
        if skipFlag
            continue
        end
        
        %Compute distance and save
        d = distance_fun(c_i,c_j);
        if d < R
            matches = matches + 1;
            motifs{i}(matches) = j;
        end
    end
    
    %Remove NaNs and include current subsequence (motif center)
    motifs{i}(isnan(motifs{i})) = [];
    if ~isempty(motifs{i});
        motifs{i} = [i motifs{i}];
    end
end

%Sort by largest motifs
[~,idx] = sort(cellfun(@length, motifs), 'descend');
motifs = motifs{idx};
%Remove empty motifs
motifs(cellfun(@isempty, motifs)) = [];

%Display results:
if displayFlag
    for i = 1:3
        subplot(3,2,2*i-1)
        title([num2str(i), '-Motif of length', num2str(l)])
        plot(timeSeries(motifs{i}(1):motifs{i}(1)+l), 'linewidth', 2);
        
        subplot(3,2,2*i)
        title(['Sequences following motif ', num2str(i)])
        for j = 1:length(motifs{i})
            hold on
            plot(timeSeries(motifs{i}(j):motifs{i}(j)+l))
        end
    end
end

end
