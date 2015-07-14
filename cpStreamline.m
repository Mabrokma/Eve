%Average multiple CompiledParticles structures and return traces for raw,
%unbinned fluo, binned mean fluo, fraction of transcribing nuclei/bin, and
%total mRNA produced, binned and unbinned.  Each strucure will have third
%dimension 3, where the field(:,:,1) are means, field(:,:,2) are std, and
%field(:,:,3) is number of embryos averaged 
%Additionally, the indices of the consensus AP region will be given as its
%own field, i.e. the region in which all embryos had data


%Averaging will come later, for now this only does one embryo at a time
function miniCP = cpStreamline(cp, varargin)
startFrame = 1;
for i = 2:nargin
    if strcmpi(varargin{i}, 'nc14')
        startFrame = cp.nc14;
    end
end

miniCP = struct;

%Size of arrays to extract, remove irrelevant particles if necessary
traceLength = length(cp.ElapsedTime) - startFrame + 1;
whichParticles = [cp.CompiledParticles.firstFrame] >= startFrame;
cp.CompiledParticles = cp.CompiledParticles(whichParticles);
nParticles = lenght(cp.CompiledParticles);

%Binning information
nBins = length(cp.APbinID)-1;
binSize = cp.APbinID(2);

%Create trace of position parallel to AllTracesVector
posTrace = NaN(traceLength,nParticles);
fluoTrace = NaN(traceLength,nParticles);

for t = startFrame:traceLength
    [posInFrame, fluoInFrame, particlesInFrame] = ...
        getParticlesInFrame(cp.CompiledParticles,t);
    
    if ~isempty(particlesInFrame)
        posTrace(t, particlesInFrame(:,1)) = posInFrame;
        fluoTrace(t, particlesInFrame(:,1)) = fluoInFrame;
    end
end

%Interpolate holes in the data
for i = 1:nParticles
    for t = 2:traceLength-1
        if isnan(fluoTrace(t,i)) ...
                && ~isnan(fluoTrace(t+1,i)) ...
                && ~isnan(fluoTrace(t-1,i))
            fluoTrace(t,i) = ...
                (fluoTrace(t+1,i) + fluoTrace(t-1,i))/2;
            posTrace(t,i) = ...
                (posTrace(t+1,i) + posTrace(t-1,i)) / 2;
        end
    end
end

%Initialize traces for nc14 & binning
activeNuclei = zeros(traceLength,nBins,3);
totalNuclei = zeros(traceLength,nBins,2);
meanFluo = NaN(traceLength,nBins,3);
%for mRNA, sum columns of fluoTrace treating NaNs as zeros, return zeros
%that persist to NaNs
mRNATrace = fluoTrace;
mRNATrace(isnan(mRNATrace)) = 0;
mRNATrace = cumsum(mRNATrace);
mRNATrace(mRNATrace == 0) = NaN;

meanmRNA = NaN(traceLength,nBins,3);

%Bin by mean AP position--make this better by allowing particles to drift
%TODO - change positions of particles everywhere to reflect positions of
%the nucleus associated with the particle
meanAP = nanmean(posTrace);

[~, ~, whichBin] = histcounts(meanAP, cp.APbinID);
miniCP.minAP = min(whichBin(whichBin~=0));
miniCP.maxAP = max(whichBin);
miniCP.meanAP = meanAP;

for x = 1:nBins
    for t = 1:traceLength
        %Get number of particles in current AP bin and time point and the
        %indices to which they correspond
        [nP, ~, which] = histcounts(posTrace(t,:),...
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
        totalNuclei(t,x,1) = nN1;
        totalNuclei(t,x,2) = nN2;
        %Get mean fluorescence for bin and time point via fluoTrace
        meanFluo(t,x,1) = nanmean(fluoTrace(t,logical(which)));
        meanFluo(t,x,2) = nanstd(fluoTrace(t,logical(which)),1);
        meanFluo(t,x,3) = nP;
    end
end




%Export to the final structure
miniCP.fluoTrace = fluoTrace;
miniCP.posTrace = posTrace;
miniCP.activeNuclei = activeNuclei;
miniCP.totalNuclei = totalNuclei;
miniCP.meanFluo = meanFluo;
miniCP.nc14 = cp.nc14;
miniCP.ElapsedTime = cp.ElapsedTime;

if isfield(cp.CompiledParticles, 'whichStripe') && startFrame > 1
    miniCP.whichStripe = [cp.CompiledParticles.whichStripe];
    miniCP.d2Centroid = [cp.CompiledParticles.d2Centroid];
end
