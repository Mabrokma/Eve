%Script to analyze all the data and generate the figures I want
clear

%Load all of the original CompiledParticles structures
GFP(3) = struct('nc14',[], 'ElapsedTime',[], 'CompiledParticles',[],...
        'AllTracesVector',[],'EllipsesFilteredPos',[],'EllipsePos',[],...
        'APbinID',[],'NEllipsesAP',[]);
GFPbyStripe = cell(1,3);
miniGFP(3) = struct('minAP', [], 'maxAP', [], 'meanAP', [],...
    'fluoTrace', [], 'posTrace', [], 'activeNuclei', [], 'totalNuclei', [], 'meanFluo', [],...
    'nc14', [], 'ElapsedTime', [], 'whichStripe', [], 'd2Centroid', []);
noGFP(3) = struct('nc14',[], 'ElapsedTime',[], 'CompiledParticles',[],...
        'AllTracesVector',[],'EllipsesFilteredPos',[],'EllipsePos',[],...
        'APbinID',[],'NEllipsesAP',[]);
noGFPbyStripe = cell(1,3);
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

h = waitbar(0);
for i = 1:3
    waitbar((6*i-5)/18,h)
    GFP(i) = load(GFPnames{i}, 'nc14', 'ElapsedTime', 'CompiledParticles',...
        'AllTracesVector','EllipsesFilteredPos','EllipsePos',...
        'APbinID','NEllipsesAP');
    waitbar((6*i-4)/18,h)

    noGFP(i) = load(noGFPnames{i}, 'nc14', 'ElapsedTime', 'CompiledParticles',...
        'AllTracesVector','EllipsesFilteredPos','EllipsePos',...
        'APbinID','NEllipsesAP');

    waitbar((6*i-3)/18,h)
    [GFPbyStripe{i}, GFP(i)] = sortByStripe(GFP(i));
    waitbar((6*i-2)/18,h)
    [noGFPbyStripe{i}, noGFP(i)] = sortByStripe(noGFP(i));
    
    waitbar((6*i-1)/18,h)
    miniGFP(i) = cpStreamline(GFP(i));
    waitbar((6*i)/18,h)
    mininoGFP(i) = cpStreamline(noGFP(i));
end
close(h)
clear GFPnames noGFPnames
