function [pos, fluo, indices] = getParticlesInFrame(CP, frame)

%Search particles of structure for current frame
pos = [];
fluo = [];
indices = [];
for i = 1:length(CP)
    thisParticle = CP(i);
    %determine if particle exists in current frame
    if ~isempty(find(thisParticle.Frame==frame,1))
        %extract index of current frame
        thisFrame = find(thisParticle.Frame==frame);
        %Save coordinates of particle and frame for reference later
        % CP(indices(:,1)) is all particles in the frame
        % CP(indices(n,1).FIELD(indices(n,2))) is the value of the field
        % for the given particle and frame
        indices(end+1,1) = i;
        indices(end,2) = thisFrame;
        
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
