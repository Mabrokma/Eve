function [genotype, avgGenotype] = loadGenotype(n)
%**************************************************************************
%This is the core format i'm using for post-analysis right now, each embryo
%gets a 'genotype' struct, similar genotypes get stored in a struct array
%for that genotype- this contains the original CompiledParticles structure
%for reference if necessary, as well as the streamlined trace vectors
%extracted by the other functions
%
%Load in a series of n CompiledParticles structures, extract the relevant
%traces, and spit out the results as a single struct as well as the average
%trace of the entire genotype.  The original compiled particles structure
%is saved in genotype.CP, and the prefix in genotype.Name.  The individual
%traces are described in extractTraces()
%
%Relies on uigetdir to choose the locations of the CompiledParticles
%structures and direct input to describe the genotype.  descriptions should
%not include spaces as they may be used in paths / as filenames
%
% Dependencies: extractTraces.m, averageTraces.m
% RW 7/2015
%**************************************************************************

genotype(n) = struct('Name',[],'desc',[],'CP',[],...
    'rawTraces',[],'standardTraces',[],'binTraces',[],'binAllTraces',[]);

for i = 1:n
    loadPath = [uigetdir, '/CompiledParticles.mat'];
    slashes = strfind(loadPath, filesep);
    genotype(i).Name = loadPath(slashes(end-1)+1:slashes(end)-1);
    genotype(i).desc = input('Input genotype descrption:\n?>','s');
    genotype(i).CP = load(loadPath);
    [genotype(i).rawTraces, genotype(i).standardTraces, ...
        genotype(i).binTraces, genotype(i).binAllTraces] = ...
        extractTraces(genotype(i).CP);
end

avgGenotype = averageTraces(genotype);