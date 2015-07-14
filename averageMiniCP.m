function avgCP = averageMiniCP(miniCP)
avgCP = struct;

l = length(miniCP);
traceLength = inf;

for i = 1:l
    if length(miniCP(i).meanFluo) < traceLength
        traceLength = length(miniCP(i).meanFluo);
    end
end


meanFluoStack = NaN(traceLength,100,l);
activeNucleiStack = NaN(traceLength,100,l);

for i = 1:l
    if i == 1
        avgCP.meanAP = miniCP(i).meanAP;
        avgCP.fluoTrace = miniCP(i).fluoTrace(1:traceLength,:);
        avgCP.posTrace = miniCP(i).posTrace(1:traceLength,:);
    else
        avgCP.meanAP = [avgCP.meanAP, miniCP(i).meanAP];
        avgCP.fluoTrace = ...
            [avgCP.fluoTrace, miniCP(i).fluoTrace(1:traceLength,:)];
        avgCP.fluoTrace = ...
            [avgCP.posTrace, miniCP(i).posTrace(1:traceLength,:)];
    end
    
    meanFluoStack(1:traceLength,:,i) = ...
        miniCP(i).meanFluo(1:traceLength,:,1);
    activeNucleiStack(1:traceLength,:,i) = ...
        miniCP(i).activeNuclei(1:traceLength,:,1);

end

avgCP.minAP = min([miniCP.minAP]);
avgCP.maxAP = max([miniCP.maxAP]);

avgCP.meanFluo(:,:,1) = nanmean(meanFluoStack(1:traceLength,:,:),3);
avgCP.activeNuclei(:,:,1) = nanmean(activeNucleiStack(1:traceLength,:,:),3);

avgCP.meanFluo(avgCP.meanFluo==0) = NaN;
avgCP.activeNuclei(avgCP.activeNuclei==0) = NaN;