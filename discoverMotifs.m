function [motifs, distances] = ...
    discoverMotifs(timeSeries, l, R, distance_fun)
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
% argument.  Distances between subsequences appearing in the same motif are
% also saved so that filtering by a smaller R can occur at a later point.
%
% TODO - don't save all motifs, only top 5 (or something), dynamically
% check whether a given motif has enough instances to be in the top 5
%
%                      -----------------------------
%
% THIS IS SO SLOW IT HURTS, NEED TO FIX and/or optimize radius somehow
% before running the full thing
%
% ANOTHER IDEA - run on small dataset, identify potential motifs, then run
% another seeded motif discovery on the larger dataset 
%**************************************************************************

if nargin < 3
    %Default metric is taxicab
    distance_fun = @(x,y)sum(abs(minus(x,y)));
    R = 2E3*l;
elseif nargin < 4
    distance_fun = @(x,y)sum(abs(minus(x,y)));
end

displayFlag = 1;

%Limits on sequence lengths
traceLength = size(timeSeries,1);
totalLength = numel(timeSeries);

emptyThreshold = 0.6;

motifs = cell(1,totalLength - l +1);
distances = cell(1,totalLength - l + 1);

h = waitbar(0);
for i = 1:totalLength - l + 1
    progress = i/(totalLength-l + 1);
    waitbar(progress, h, ['Finding motifs...    ', ...
            num2str(progress*100,'%4.2f'), '% Complete']);
        
    %Initialize count & preallocate index matrix
    matches = 0;
    motifs{i} = NaN(1,floor(totalLength/l - traceLength*l));
    distances{i} = NaN(1,floor(totalLength/l - traceLength*l));
    
    %Subsequence i(l)
    iIdx = i:i+l-1;            %indices of subsequence i
    c_i = timeSeries(iIdx);    %values of subsequence i
    
    %Skip bad subsequences
    skipFlag = (any(mod(iIdx,traceLength) == 0) ...
        && any(mod(iIdx,traceLength) == 1)) ...  cross traces
        || sum(c_i == 0)/l > emptyThreshold;%    empty
    
    if skipFlag
        motifs{i} = [];
        continue
    end
    
    for j = 1:totalLength - l + 1;
        
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
            %Only allow one overlapping instance per motif (the best)
            if any(j == motifs{i}(end):motifs{i}(end)+l-1)
                %Replace multiple overlapping matches with the best match
                if d < distances{i}(end)
                    motifs{i}(matches) = j;
                    distances{i}(matches) = d;
                else
                    continue
                end
            else
                %Nonoverlapping matches are appended
                matches = matches + 1;
                motifs{i}(matches) = j;
                distances{i}(matches) = d;
            end
        end
    end
    
    %Remove NaNs and include current subsequence (motif center)
    motifs{i}(isnan(motifs{i})) = [];
    distances{i}(isnan(distances{i})) = [];
    if ~isempty(motifs{i});
        motifs{i} = [i motifs{i}];
        distances{i} = [0 distances{i}];
    end
end

waitbar(0.99, h, 'Sorting...')
%TODO - filter out overlapping motifs, i.e. d(motif1, motif2) > 2R

%Sort by largest motifs
[~,idx] = sort(cellfun(@length, motifs), 'descend');
motifs = motifs(idx);
distances = distances(idx);
%Remove empty motifs
motifs(cellfun(@isempty, motifs)) = [];
distances(cellfun(@isempty, distances)) = [];

delete(h);
%Display results:
if displayFlag
    for i = 1:3
        %Plot all instances of a motif together
        subplot(3,2,2*i-1)
        hold on
        plot(timeSeries(motifs{i}(1):motifs{i}(1)+l), 'linewidth', 2);
        for j = 2:length(motifs{i})
            hold on
            plot(timeSeries(motifs{i}(j):motifs{i}(j)+l))
        end
        title([num2str(i), '-Motif of length', num2str(l)])
        xlabel('Relative time')
        ylabel('Fluorescence')
        
        %Plot instances of motif as they appear
        subplot(3,2,2*i)
        for j = 1:length(motifs{i})
            hold on
            plot(mod(motifs{i}(j):(motifs{i}(j)+l),traceLength),...
                timeSeries(motifs{i}(j):motifs{i}(j)+l))
        end
        title('Real-time occuence of motif')
        xlabel('Absolute time')
        ylabel('Fluorescence')
    end
end

end
