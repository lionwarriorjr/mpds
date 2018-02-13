function [header, feature] = extractStreamBattV1(data)

[m n] = size(data);
sum_sec = sum(data(:,1));
batt = data(1,2) - data(m,2);

header = {'interval_sec', 'batt_level', 'batt_level_per_sec'};
feature = [sum_sec batt batt/sum_sec];

