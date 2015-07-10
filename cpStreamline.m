%Average multiple CompiledParticles structures and return traces for raw,
%unbinned fluo, binned mean fluo, fraction of transcribing nuclei/bin, and
%total mRNA produced, binned and unbinned.  Each strucure will have third
%dimension 3, where the field(:,:,1) are means, field(:,:,2) are std, and
%field(:,:,3) is number of embryos averaged 
%Additionally, the indices of the consensus AP region will be given as its
%own field, i.e. the region in which all embryos had data
%this will also suck out only nc14


%Averaging will come later, for now this only does one embryo at a time
function shortCP = cpStreamline(cp)
shortCP = struct;

%Size of arrays to extract
traceLength = length(cp.ElapsedTime);
nParticles = length(cp.CompiledParticles);
%Binning information
nBins = length(cp.APbinID)-1;
binSize = cp.APbinID(2);

%Create trace of position parallel to AllTracesVector
posTrace = NaN(traceLength,nParticles);

for t = 1:traceLength
    [posInFrame, ~, particlesInFrame] = ...
        getParticlesInFrame(cp.CompiledParticles,t);
    
    if ~isempty(particlesInFrame)
        posTrace(t, particlesInFrame(:,1)) = posInFrame;
    end
end

%Extract relevant subset of particles
for i = 1:length(cp.CompiledParticles)
    firstParticle = i;
    
    if cp.CompiledParticles(i).FirstFrame > cp.nc14
        break
    end
end

%Indices for particles/time of interest
time = cp.nc14:length(cp.ElapsedTime);
particles = firstParticle:length(cp.CompiledParticles);
shortTraceLength = length(time);
shortnParticles = length(particles);

shortFluoTrace = cp.AllTracesVector(time,particles);
shortPosTrace = posTrace(time,particles);

%Interpolate holes in the data
for i = 1:shortnParticles
    for t = 2:shortTraceLength-1
        if isnan(shortFluoTrace(t,i)) ...
                && ~isnan(shortFluoTrace(t+1,i)) ...
                && ~isnan(shortFluoTrace(t-1,i))
            shortFluoTrace(t,i) = ...
                (shortFluoTrace(t+1,i) + shortFluoTrace(t-1,i))/2;
            shortPosTrace(t,i) = ...
                (shortPosTrace(t+1,i) + shortPosTrace(t-1,i)) / 2;
        end
    end
end

%Initialize traces for nc14 & binning
activeNuclei = zeros(shortTraceLength,nBins);
meanFluo = NaN(shortTraceLength,nBins,3);
mRNATrace = NaN(shortTraceLength,shortnParticles,3);
meanmRNA = NaN(shortTraceLength,nBins,3);

%Bin by mean AP position--make this better by allowing particles to drift
meanAP = nanmean(shortPosTrace);

[~, ~, whichBin] = histcounts(meanAP, cp.APbinID);
shortCP.minAP = min(whichBin(whichBin~=0))-1;
shortCP.maxAP = max(whichBin)+1;
shortCP.meanAP = meanAP;

for x = 1:nBins
    for t = 1:shortTraceLength
        %Get number of particles in current AP bin and time point and the
        %indices to which they correspond
        [nP, ~, which] = histcounts(shortPosTrace(t,:),...
            [cp.APbinID(x), cp.APbinID(x+1)]);
        %Get number of (approved) nuclei in AP bin and time point
        nN1 = histcounts(cp.EllipsesFilteredPos{t},...
            [cp.APbinID(x), cp.APbinID(x+1)]);
        nN2 = histcounts(cp.EllipsePos{t},...
            [cp.APbinID(x), cp.APbinID(x+1)]);
        %Calculate fraction of active nuclei in bin and time point
        activeNuclei(t,x,1) = nP/nN1;
        activeNuclei(t,x,2) = nP/nN2;
        activeNuclei(t,x,3) = nP;
        %Get mean fluorescence for bin and time point via fluoTrace
        meanFluo(t,x,1) = nanmean(shortFluoTrace(t,logical(which)));
        meanFluo(t,x,2) = nanstd(shortFluoTrace(t,logical(which)),1);
        meanFluo(t,x,3) = nP;
    end
end


%Export to the final structure
shortCP.fluoTrace = shortFluoTrace;
shortCP.posTrace = shortPosTrace;
shortCP.activeNuclei = activeNuclei;
shortCP.meanFluo = meanFluo;

if isfield(cp.CompiledParticles, 'whichStripe')
    shortCP.whichStripe = [cp.CompiledParticles(particles).whichStripe];
    shortCP.d2Centroid = [cp.CompiledParticles(particles).d2Centroid];
end
