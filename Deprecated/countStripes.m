function [nStripes, foundStripes] = countStripes(CP, syncedFrames)
%Expand to determine transitionFrame and startFrame

for i = 10:length(syncedFrames);
    %Smooth and find peaks, determine how many stripes there are by the
    %number of peaks in the smoothed histogram
    [APpos, ~] = getParticlesInFrame(CP, syncedFrames(i));
    posBins = 0:0.01:1;
    posCounts = histcounts(APpos, posBins);
    w = 3; %moving window size
    smoothCounts = conv(posCounts, ones(1,w)/w, 'same');
    foundStripes(i) = length(findpeaks(smoothCounts));
end

nStripes = mode(foundStripes);