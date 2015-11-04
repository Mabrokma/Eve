function [] = plotStripes(traces, stripes)

close all;
figure('Units', 'Normalized', 'Position', [0 0 1 1]);
edges = 0.2:0.01:0.85;
for t = 1:size(traces,1)
    %Histogram of entire pattern
    subplot(2,1,1)
    cla; hold on;
    [counts] = histcounts(traces(t,:,3),edges);
    plot(edges(1:end-1)+0.005, counts);
    axis([0.2 0.85 0 10])
    
    %Gaussian fit for each stripe
    subplot(2,1,2)
    cla; hold on;
    for s = 1:max(stripes)
        try
            pd = fitdist(traces(t,stripes==s,3)',...
                'Normal');
            pos = min(edges):0.002:max(edges);
            freq = pdf(pd, pos);
            plot(pos,freq);
            axis([0.2 0.85 0 40])
        catch
            %Empty ones
        end
    end
    pause(0.2)
end