%Load relevant CompiledParticles data structures

function [E, C, syncedFrames, nParticles, nFrames] = ...
    loadCP(traceLength, varargin)

%Parse args
if isempty(varargin)
    n = 1;
else
    n = varargin{1};
end

%Load relevant CompiledParticles data strucstures
for k = 1:n
    loadPath = uigetdir(cd, ['Choose processed data set']);
    try
        E{k} = load([loadPath,filesep, 'CompiledParticles.mat'],...
            'CompiledParticles', 'ElapsedTime', 'nc14');
    catch
        error('CompiledParticles.mat not found')
    end
    
    %Get name of dataset
    dashes=strfind(loadPath,filesep);
    Prefix{k}=loadPath((dashes(end)+1):end);
    %Shorten name for convenience
    C{k} = E{k}.CompiledParticles;
    %Get number of particles and length of each trace
    nParticles(k) = length(C{k});
    nFrames(k) = length(E{k}.ElapsedTime);
    
    %Synchronize start point as beginning of nc14, take [movieLength] frames
    syncedFrames(k,:) = (E{k}.nc14):(E{k}.nc14+traceLength-1);
    
    %TODO make number of frames customizable and/or fixed to length of
    %input traces.  also synchronize based on time and not frames
end