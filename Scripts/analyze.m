%Script to analyze all the data and generate the figures I want
tic; parpool;
clear GFP noGFP miniGFP mininoGFP;

%Load all of the original CompiledParticles structures
GFP(3) = struct('nc14',[], 'ElapsedTime',[], 'CompiledParticles',[],...
        'AllTracesVector',[],'EllipsesFilteredPos',[],'EllipsePos',[],...
        'APbinID',[]);
GFPbyStripe = cell(1,3);
miniGFPbyStripe = cell(1,3);
miniGFP(3) = struct('minAP', [], 'maxAP', [], 'meanAP', [],...
    'fluoTrace', [], 'posTrace', [], 'activeNuclei', [], 'totalNuclei', [], 'meanFluo', [],...
    'nc14', [], 'ElapsedTime', [], 'whichStripe', [], 'd2Centroid', []);
noGFP(3) = struct('nc14',[], 'ElapsedTime',[], 'CompiledParticles',[],...
        'AllTracesVector',[],'EllipsesFilteredPos',[],'EllipsePos',[],...
        'APbinID',[]);
noGFPbyStripe = cell(1,3);
mininoGFPbyStripe = cell(1,3);
mininoGFP(3) = struct('minAP', [], 'maxAP', [], 'meanAP', [],...
    'fluoTrace', [], 'posTrace', [], 'activeNuclei', [], 'totalNuclei', [], 'meanFluo', [],...
    'nc14', [], 'ElapsedTime', [], 'whichStripe', [], 'd2Centroid', []);

%load them all
filePath = '/Users/onetrob/Documents/Princeton/Gregor/LivemRNA/Data/DynamicsResults/2015-03-';
GFPnames = {[filePath, '10-eveGFP_Embryo14/', 'CompiledParticles.mat'],...
    [filePath, '11-eveGFP_Embryo02/', 'CompiledParticles.mat'],...
    [filePath, '12-eveGFP_Embryo05/', 'CompiledParticles.mat']};
noGFPnames = {[filePath, '11-eveNoGFP_Embryo03/', 'CompiledParticles.mat'],...
    [filePath, '12-eveNoGFP_Embryo03/', 'CompiledParticles.mat'],...
    [filePath, '13-eveNoGFP_Embryo04/', 'CompiledParticles.mat']};

h = waitbar(0, 'Total Progress');
parfor i = 1:3
    GFP(i) = load(GFPnames{i}, 'nc14', 'ElapsedTime', 'CompiledParticles',...
        'AllTracesVector','EllipsesFilteredPos','EllipsePos',...
        'APbinID');

    noGFP(i) = load(noGFPnames{i}, 'nc14', 'ElapsedTime', 'CompiledParticles',...
        'AllTracesVector','EllipsesFilteredPos','EllipsePos',...
        'APbinID');

    waitbar((10*i-9)/30,h)
    [GFPbyStripe{i}, GFP(i)] = sortByStripe(GFP(i));
    [noGFPbyStripe{i}, noGFP(i)] = sortByStripe(noGFP(i));
    
    waitbar((10*i-8)/30,h)
    miniGFP(i) = cpStreamline(GFP(i), 'nc14');
    waitbar((10*i-7)/30,h)
    mininoGFP(i) = cpStreamline(noGFP(i), 'nc14');
    
    %Do more 
    miniGFPbyStripe{i} = cell(1,7);
    mininoGFPbyStripe{i} = cell(1,7);
    for s = 1:7
        waitbar((10*i-7+s)/30,h)
        try
            miniGFPbyStripe{i}{s} = ...
                cpStreamline(GFPbyStripe{i}{s},'nc14');
            mininoGFPbyStripe{i}{s} = ...
                cpStreamline(noGFPbyStripe{i}{s},'nc14');
        catch ME
            fprintf('Skipping Stripe %d, no data\n', s)
            warning(ME.message)
            for l = 1:length(ME.stack)
                display(['In ', ME.stack(l).name, ' Line: ',...
                    num2str(ME.stack(l).line)])
            end
        end
    end
end
close(h)
clear GFPnames noGFPnames i s h ME

save('Processed.mat')
time = toc;
fprintf('\n\nTime Elapsed: %3.1f minutes\n', time/60);