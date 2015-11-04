function [genotype] = loadGenotype(n)
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

%If number of sets isn't specified, allow user to search for them
if nargin == 0
    %Data must fall in one folder (for now)
    searchPath = uigetdir('Please select folder containing data');
    datasets = dir(searchPath);
    fnames = {datasets.name};
    nTot = length(fnames);
    
    %Get regexp for searching
    query = input('Specify filenames via regular expresion:\n?>','s');
    matches = regexp(fnames,query);
    %Filter
    fnames = fnames(~cellfun(@isempty,matches));
    n = length(fnames);
    
    fprintf('\n** %i/%i Datasets match search **\n\n', n, nTot)
    clear matches query datasets
end

%Initialize Data structure
genotype(n) = struct('Name',[],'desc',[],'CP',[],...
    'rawTraces',[],'standardTraces',[],'binTraces',[],'binAllTraces',[]);

desc = input('Input genotype descrption:\n?>','s');

for i = 1:n
    %Locate CompiledParticles
    if exist('fnames', 'var')
        loadPath = [searchPath, filesep, fnames{i},...
            filesep, 'CompiledParticles.mat'];
    else
        %Get each folder individually if unspecified
        loadPath = [uigetdir, filesep, 'CompiledParticles.mat'];
    end
    
    %Attempt to load CompiledParticles
    try
        genotype(i).CP = load(loadPath);
    catch
        %Skip bad ones
        warning(['No CompiledParticles found in dataset #', num2str(i)])
        continue
    end
    
    %Get prefix name and description
    slashes = strfind(loadPath, filesep);
    genotype(i).Name = loadPath(slashes(end-1)+1:slashes(end)-1);
    genotype(i).desc = desc;
    
    %Make trace vectors
    [genotype(i).rawTraces, genotype(i).binTraces] = ...
        extractTraces(genotype(i).CP);
    [genotype(i).nc14, genotype(i).bin14,...
        genotype(i).nc13, genotype(i).bin13] = ...
        standardizeTraces(genotype(i));
end