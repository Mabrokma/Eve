%Sort CompiledParticles structures into seven identical structures that
%each represent a single eve stripe.  These can be fed into any
%post-processing code as if it were a true CompiledParticles structure
%
%

%The arg cp is a CompiledParticles.mat file
function [cpByStripe, centroids] = sortByStripe(cp)

DEFAULT_CENTROIDS = [0.3 0.38 0.46 0.525 0.59 0.67 0.76]';
nParticles = length(cp.CompiledParticles);
%Each entry in the cell will be a CompiledParticles-style structure for the
%corresponding eve stripe
cpByStripe = cell(3,7);

%Filter out only nc14
firstFrame = cp.nc14;
frameShift = firstFrame - 1;
lastFrame = length(cp.ElapsedTime);
frameOfnc14 = firstFrame:lastFrame;
nFrames = length(frameOfnc14);

%Determine which stripes are present in the viewing window
%TODO - make these frames determined dynamically
earlyRange = 1:45;
clusterRange = 46:85;  %use kmeans in this range
lateRange = 86:(lastFrame - frameShift);


%Maybe make this bit a separate function?
potentialCentroids = cell(1, length(clusterRange));
foundStripes = NaN(length(clusterRange),1);
left = NaN(length(clusterRange),1);
APpos = cell(1,length(clusterRange));
binSize = cp.APbinID(2);

%NOTE: APpos is indexed from 1:length(clusterRange) but the values of
%clusterRange do not begin at 1; t corresponds to frame, T to index
T = 0;
for t = clusterRange
    %Use capital T for indexing 
    T = T+1;
    APpos{T} = getParticlesInFrame(cp.CompiledParticles, frameOfnc14(t));
    
    %Make histogram of frame, smooth, count peaks as proxy for stripes
    posCounts = histcounts(APpos{T}, cp.APbinID);
    w = 4; %moving window size
    smoothCounts = conv(posCounts, ones(1,w)/w, 'same'); %filter
    [~,potentialCentroids{T}] = findpeaks(smoothCounts);
    potentialCentroids{T} = potentialCentroids{T}*binSize - binSize/2;
    foundStripes(t) = length(potentialCentroids{T});
    
    distancesFromCentroids = ...
        abs(potentialCentroids{T}(1)-DEFAULT_CENTROIDS);
    %Compare found centroids to default centroids
    [d, left(T)] = min(distancesFromCentroids);
end
nStripes = mode(foundStripes);
firstStripe = mode(left);
lastStripe = firstStripe+nStripes-1;
stripeShift = firstStripe - 1;

fprintf('Stripes detected: %d\n', nStripes);
if lastStripe > 7
    fprintf('Trimming extraneous stripes\n');
    nStripes = length(firstStripe:7);
    lastStripe = 7;
end

stripeNumber = firstStripe:lastStripe;

%First use CompiledParticles.APpos to perform clustering on particles alone
%Cluster in the middle of nc14, only with particles that are active at that
%time, then do assignment based on mean AP positions
%i.e. each particle gets assigned to only one cp structure
fprintf('Clustering on stripes %d - %d\n', firstStripe, lastStripe);
centroids = NaN(nFrames, nStripes);
T = 0;
for t = clusterRange
    %Captial T for indexing w/i clusterRange
    T = T+1;
    %Do clustering in the mid-nc14 period while stripes are well-defined
    [~, centroids(t,:)] = kmeans(APpos{T}', nStripes, ...
        'Start', DEFAULT_CENTROIDS(firstStripe:lastStripe));
    centroids(t,:) = sort(centroids(t,:));
end
%Set the centroids for the edges
earlyCentroids = mean(centroids(clusterRange(1:3),:));
centroids(earlyRange,:) = repmat(earlyCentroids,length(earlyRange),1);

lateCentroids = mean(centroids(clusterRange(end-2:end),:));
centroids(lateRange,:) = repmat(lateCentroids,length(lateRange),1);

%Curate the centroids here
meanCentroids = mean(centroids);
fprintf('\nStripe Number | Default Centroid | Mean Centroid\n')

for s = stripeNumber
    fprintf('      %1d       |       %3.2f       |      %3.2f     \n',...
        s, DEFAULT_CENTROIDS(s), meanCentroids(s-firstStripe+1))
end

%Filter centroids that fall larger than a certain tolerance from the mean
%centroid of that stripe.  Replace it with the average of the two adjacent
tolerance = 0.02;
for s = 1:nStripes
    %edge
    t = clusterRange(1);
    if centroids(t,s) - meanCentroids(s) > tolerance
        bad = centroids(t,s);
        centroids(t,s) = centroids(t+1,s);
        fprintf('Intolerable centroid in Frame: %d, Stripe: %d\n',...
            t,stripeNumber(s))
        fprintf('Default: %3.2f, Mean: %3.2f, Old: %3.2f, New: %3.2f\n',...
            DEFAULT_CENTROIDS(stripeNumber(s)), meanCentroids(s),...
            bad, centroids(t,s))
    end
    %center
    for t = clusterRange(2:end-1)
        if abs(centroids(t,s) - meanCentroids(s)) > tolerance
            bad = centroids(t,s);
            centroids(t,s) = (centroids(t+1,s) + centroids(t-1,s))/2;
            fprintf('Intolerable centroid in Frame: %d, Stripe: %d\n',...
                t,stripeNumber(s))
            fprintf('--> Default: %3.2f, Mean: %3.2f, Old: %3.2f, New: %3.2f\n',...
                DEFAULT_CENTROIDS(stripeNumber(s)), meanCentroids(s),...
                bad, centroids(t,s))
        end
    end
    %edge
    t = clusterRange(end);
    if centroids(t,s) - meanCentroids(s) > tolerance
        bad = centroids(t,s);
        centroids(t,s) = centroids(t-1,s);
        fprintf('Intolerable centroid in Frame: %d, Stripe: %d\n',...
            t,stripeNumber(s))
        fprintf('Default: %3.2f, Mean: %3.2f, Old: %3.2f, New: %3.2f\n',...
            DEFAULT_CENTROIDS(stripeNumber(s)), meanCentroids(s),...
            bad, centroids(t,s))
    end
end

%Centroids will be a vector with a value at each time point, particles will
%be assigned to stripes based on the centroid they are closest to in the
%first frame in which they appear
whichStripe = NaN(1,nParticles);
for i = 1:nParticles
    %Figure out initial frame and position of the particle
    frame = cp.CompiledParticles(i).FirstFrame;
    pos = cp.CompiledParticles(i).APpos(1);
    
    [cp.CompiledParticles(i).D2Centroid, whichStripe(i)] = ...
        min(abs(centroids(frame,:) - pos));
end

%Shift over if any stripes are missing
whichStripe = whichStripe + stripeShift;

%Assign elements of CompiledParticles based on clustering to the cell
for s = firstStripe:lastStripe
    cpByStripe{1,s} = cp.CompiledParticles(whichStripe==s);
end

%assign ellipses w/o particles as well
%initialize structure
for s = firstStripe:lastStripe
    cpByStripe{2,s} = cell(1,nFrames);
    cpByStripe{3,s} = cell(1,nFrames);
end
%Assign filtered ellipses to the second row of the cell
for t = firstFrame:lastFrame
    pos = cp.EllipsesFilteredPos{t};
    whichStripe = NaN(1,length(pos));
    %Calculate nearest stripe
    for i = 1:length(pos)
        [~, whichStripe(i)] = min(abs(centroids(t-frameShift,:) - pos(i)));
    end
    whichStripe = whichStripe + stripeShift;
    
    %Put results into structure
    for s = firstStripe:lastStripe
        cpByStripe{2,s}{t - frameShift} = pos(whichStripe==s);
    end
end

%Repeat for unfiltered ellipses
for t = firstFrame:lastFrame
    pos = cp.EllipsePos{t};
    whichStripe = NaN(1,length(pos));
    %Calculate nearest stripe
    for i = 1:length(pos)
        [~, whichStripe(i)] = min(abs(centroids(t-frameShift,:) - pos(i)));
    end
    whichStripe = whichStripe + stripeShift;
    
    %Put results into structure
    for s = firstStripe:lastStripe
        cpByStripe{3,s}{t - frameShift} = pos(whichStripe==s);
    end
end


