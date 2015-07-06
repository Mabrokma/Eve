function [meanFluoFig, activeNucleiFig] = ...
    plotAllStripes(stripes, Embryo, Prefix)

%Check args, throw errors
if ~isfield(Embryo, 'activeNuclei') || ~isfield(Embryo, 'meanFluo')
    error('Please feed me real embryos, you are a dufus')
end

traceLength = length(Embryo.activeNuclei);

%For each embryo, plot all stripes on the same axis in different colors
%Initialize figures, take up entire screen
meanFluoFig = ...
    figure('units','normalized','outerposition',[0 0 1 1]);
activeNucleiFig = ...
    figure('units','normalized','outerposition',[0 0 1 1]);

leg = cell(1,7); %Legend strings
for s = 1:7
    figure(meanFluoFig)
    hold on
    if ~isempty(Embryo.meanFluo(:,s))
        plot(1:traceLength, Embryo.meanFluo(:,s))
    else
        plot(0,0, '--')
    end
    xlabel('Frame Number')
    ylabel('Mean Fluorescence')
    title(['Fluorescence over time: ', Prefix])
    
    figure(activeNucleiFig)
    hold on
    if ~isempty(Embryo.meanFluo(:,s))
        plot(1:traceLength, Embryo.activeNuclei(:,s),...
            'linewidth', 2)
    else
        plot(0,0, '--')
    end
    xlabel('Frame Number')
    ylabel('Fraction of Transcribing Nuclei')
    title(['Transcribing nuclei over time: ', Prefix])
    leg{s} = ['Stripe #', num2str(s)];
end

%Set legends
figure(meanFluoFig)
legend(leg, 'Location', 'Best')
figure(activeNucleiFig)
legend(leg, 'Location', 'Best')