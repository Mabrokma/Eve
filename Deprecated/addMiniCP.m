%Add two miniCP structures together for the purpose of modeling an enhancer
%that controls more than one stripe

function newCP = addMiniCP(cp1,cp2)
newCP = struct;
newCP.minAP = min(cp1.minAP,cp2.minAP);
newCP.maxAP = max(cp1.maxAP,cp2.maxAP);

newCP.activeNuclei = cp1.activeNuclei + cp2.activeNuclei;
newCP.ElapsedTime = cp1.ElapsedTime;
newCP.totalNuclei = cp1.totalNuclei;