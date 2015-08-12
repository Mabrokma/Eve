function [whichParticle, whichFrame, varargout] = ...
    particlesInFrame(CP, t, varargin)
%**************************************************************************
%Find out which particles exist in a given time point, t, and return the
%proper indices for accessing these particles in the CompiledParticles
%Structure.
%
%Also return the corresponding values fields of CompiledParticles supplied
%as additional args in the form of string names of fields.
%
%CP(particles) is all particles in the frame
%CP(particles(t).FIELD(frame(t))) is the value of the field for the given 
%particle and frame
%
% Dependencies: none
% RW 7/2015
%**************************************************************************

%Preallocate the maximum possible size
whichParticle = NaN(1,length(CP));
whichFrame = NaN(1,length(CP));
n = 1; %Counter of particles in frame

%Check that requested fields exist
if nargin>2
    nFields = length(varargin);
    varargout = cell(1,nFields);
    for j = 1:nFields
        if isfield(CP, varargin{j})
            varargout{j} = NaN(1,length(CP));
        else
            error(['Field ', varargin{j}, ' not found'])
        end
    end
end

%Search particles of structure for current frame
for i = 1:length(CP)
    thisParticle = CP(i);
    %determine if particle exists in current frame
    if any(thisParticle.Frame == t)
        %extract index of current frame
        thisFrame = find(thisParticle.Frame==t);
        
        %Save coordinates of particle and frame
        whichParticle(n) = i;
        whichFrame(n) = thisFrame;
        n = n+1;
        
        %Get any requested properties at this time point
        for j = 1:length(varargin)
            field = thisParticle.(varargin{j});
            if length(field) > 1
                varargout{j}(n) = field(thisFrame);
            else
                varargout{j}(n) = field;
            end
        end
    end
end

%Trim extra space
whichParticle(isnan(whichParticle)) = [];
whichFrame(isnan(whichFrame)) = [];

for j = 1:length(varargin)
    varargout{j}(isnan(varargout{j})) = [];
end

end