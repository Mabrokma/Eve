function mRNA = integrateWithDegradation(binFluo, t50)
%Do a cumulative sum on the fluorescence, allowing for degradation with a
%specified half-life using first-order kinetics

binFluo(isnan(binFluo)) = 0;

if nargin < 2
    t50 = 6;
end

DEGRADATION = log(2)/t50;

mRNA = zeros(size(binFluo,1), size(binFluo,2));
mRNA(1,:) = binFluo(1,:);

for t = 2:size(binFluo,1)
    mRNA(t,:) = mRNA(t-1,:)*(1-DEGRADATION) + binFluo(t,:);
end

