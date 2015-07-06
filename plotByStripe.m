%Plot multiple traces corresponding to a single stripe
function stripeFig = ...
    plotByStripe(thisStripe, stripes, Embryos, Prefix, varargin)
%Todo - check args and throw errors
if ~isfield(Embryos, 'activeNuclei') || ~isfield(Embryos, 'meanFluo')
    error('Please feed me real embryos, you are a dufus')
end

nEmbryos = length(Embryos);
traceLength = length(Embryos(1).activeNuclei);

%Initialize individual stripe figures, take up whole screen
stripeFig = figure('units','normalized','outerposition',[0 0 1 0.5]);
%right plot - fluorescence
subplot(1,2,2)
hold on
h = zeros(1,nEmbryos);
for k = 1:nEmbryos
    if ~isempty(Embryos(k).meanFluo(:,thisStripe))
        h(k) = plot(1:traceLength, Embryos(k).meanFluo(:,thisStripe),...
            'linewidth', 2);
        
        if isfield(Embryos, 'errFluo')
            %Keep colors consistent when plotting error bars
            ax = gca;
            ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
            errorbar(1:traceLength, Embryos(k).meanFluo(:,thisStripe),...
                Embryos(k).errFluo(:,thisStripe))
        end
    else
        %Throw blank plots to keep colors consistent
        plot(0,0, '--')
    end
end
title(['Stripe #', num2str(thisStripe),...
    ': Mean fluorescence over time'])
xlabel('Frame Number')
ylabel('Mean Fluorescence')
legend(Prefix, 'Location', 'Best')

%left plot - number of transcribing nuclei
subplot(1,2,1)
hold on
for k = 1:nEmbryos
    if ~isempty(Embryos(k).activeNuclei(:,thisStripe))
        h(k) = plot(1:traceLength, Embryos(k).activeNuclei(:,thisStripe),...
            'linewidth', 2);
        if isfield(Embryos, 'errNuclei')
            %Keep colors consistent when plotting error bars
            ax = gca;
            ax.ColorOrderIndex = ax.ColorOrderIndex - 1;
            errorbar(1:traceLength, ...
                Embryos(k).activeNuclei(:,thisStripe), ...
                Embryos(k).errNuclei(:,thisStripe))
        end
    else
        %Throw blank plots to keep colors consistent
        plot(0,0, '--')
    end
end
title(['Stripe #', num2str(thisStripe),...
    ': Transcribing nuclei over time'])
xlabel('Frame Number')
ylabel('Fraction of Transcribing Nuclei')
axis([0, traceLength, 0, 1])
legend(Prefix, 'Location', 'Best')
