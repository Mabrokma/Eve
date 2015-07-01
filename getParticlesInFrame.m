function [pos, fluo] = getParticlesInFrame(CP, frame)

%Search particles of structure for current frame
pos = [];
fluo = [];
for i = 1:length(CP)
    thisParticle = CP(i);
    %determine if particle exists in current frame
    if ~isempty(find(thisParticle.Frame==frame))
        %extract index of current frame
        thisFrame = find(thisParticle.Frame==frame);
        
        %use index to find fluorescence and AP position
        thisPos = thisParticle.APposParticle(thisFrame);
        thisFluo = thisParticle.Fluo(thisFrame);
        
        %add to running structure that holds everything
        %yes this changes size every loop iteration, not sure how to
        %preallocate since final size is unknown
        pos(end+1) = thisPos;
        fluo(end+1) = thisFluo;
    end
end