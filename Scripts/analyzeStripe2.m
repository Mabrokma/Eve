%Script to analyze all the data and generate the figures I want
tic; clear Eve;

%Load all of the original CompiledParticles structures
Eve(3) = struct('nc14',[], 'ElapsedTime',[], 'CompiledParticles',[],...
        'AllTracesVector',[],'EllipsePos',[],'APbinID',[]);

miniEve(3) = struct('minAP', [], 'maxAP', [], 'meanAP', [],...
    'fluoTrace', [], 'posTrace', [], 'activeNuclei', [], 'totalNuclei', [], 'meanFluo', [],...
    'nc14', [], 'ElapsedTime', []);

%load them all
filePath = '/Users/onetrob/Documents/Princeton/Gregor/LivemRNA/Data/DynamicsResults/2014-03-';
fileNames = {[filePath, '20-Eve2E/', 'CompiledParticles.mat'],...
    [filePath, '14-Eve2B/', 'CompiledParticles.mat'],...
    [filePath, '19-Eve2A/', 'CompiledParticles.mat'],...
    [filePath, '20-Eve2B/', 'CompiledParticles.mat'],...
    [filePath, '20-Eve2C/', 'CompiledParticles.mat'],...
    [filePath, '20-Eve2D/', 'CompiledParticles.mat']};


h = waitbar(0, 'Total Progress');
h.Units = 'Normalized';
h.Position(2) = h.Position(2) + h.Position(4);
for i = 1:6
    Eve(i) = load(fileNames{i}, 'nc14', 'ElapsedTime', 'CompiledParticles',...
        'AllTracesVector','EllipsePos','APbinID');

    miniEve(i) = cpStreamline(Eve(i), 'nc14');
    waitbar(i/6,h)    
end
close(h)
clear i h filePath fileNames

save('EveStripe2.mat')
time = toc;
fprintf('\n\nTime Elapsed: %3.1f minutes\n', time/60);