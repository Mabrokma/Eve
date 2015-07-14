%Average together fluorescence of many time
function Gap = gapGeneTimeTrace(Gap, timeStep)

%Find avg AP value for each time range
timeBins = 0:timeStep:60;
oFilter = logical(Gap.Orientation == 1);
Data = Gap.Data(oFilter,:,:);
Time = Gap.Time(oFilter,:,:);

avgExp = NaN(length(timeBins)-1,1000,4);
noise = NaN(length(timeBins)-1,1000,4);

for t = 1:length(timeBins)-1
    tFilter = logical((Time > timeBins(t)) .* (Time < timeBins(t+1)));
    avgExp(t,:,:) = nanmean(Data(tFilter,:,:),1);
    noise(t,:,:) = nanvar(Data(tFilter,:,:),1,1)./avgExp(t,:,:);
    
    for i = 1:100
        Gap.avgExp(t,i,:) = nanmean(avgExp(t,(i-1)*10+1:(i*10),:),2);
    end
end


for i = 1:4
figure(i);
%subplot(1,2,1);
imagesc([0,1],[0,60], Gap.avgExp(:,:,i));
title([Gap.Name{i}, ': Expression Pattern'])
xlabel('AP Position (%EL)')
ylabel('Time into nc14 (min)')
% subplot(1,2,2);
% imagesc([0,1],[0,60], Gap.noise(:,:,i));
% title([Gap.Name{i}, ': Expression Noise'])
% xlabel('AP Position (%EL)')
% ylabel('Time into nc14 (min)')
colormap jet; colorbar;
end


%plot sum of all gap genes together?
total = Gap.avgExp(:,:,1) / max(max(Gap.avgExp(:,:,1))) ...
    + Gap.avgExp(:,:,2) / max(max(Gap.avgExp(:,:,2))) ...
    + Gap.avgExp(:,:,3) / max(max(Gap.avgExp(:,:,3))) ...
    + Gap.avgExp(:,:,4) / max(max(Gap.avgExp(:,:,4)));

figure(5)
imagesc([0,1],[0,60], total);
colormap jet; colorbar;

%Save timestep info to the returned structure
Gap.timeStep = timeStep;