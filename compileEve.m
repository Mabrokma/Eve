function eve = compileEve(Prefix)
%**************************************************************************
% Run all important functions on the stuff
%**************************************************************************
if nargin == 0
    loadPath = [uigetdir, 'CompiledParticles.mat'];
    slashes = strfind(loadPath, filesep);
    Prefix = loadPath(slashes(end-1)+1:slashes(end)-1);
else
    loadPath = ['/Users/onetrob/Documents/Princeton',...
        '/Gregor/LivemRNA/Data/DynamicsResults/', 'CompiledParticles.mat'];
end

cp = load(loadPath);

%Build structure
%This order is important so that structures may be combined
%Name of files
eve.Name = Prefix;
%One word to describe genotype
eve.desc = input('Enter description for genotype\n?>>','s');
%Original structure
eve.CP = cp;
%Raw Traces
[eve.rawTraces, eve.binTraces] = extractTraces(cp);
%Standardtraces
[eve.nc14, eve.bin14, eve.nc13, eve.bin14] = standardizeTraces(eve);