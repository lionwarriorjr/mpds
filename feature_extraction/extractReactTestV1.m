function [header, feature] = extractReactTestV1( data )
%EXTRACTACCELTESTV1 Summary of this function goes here
%   Detailed explanation goes here
[nr, nc] = size(data);
tsp = data(:,1);
btn_visiable = data(:,4);
btn_pressed = data(:,5);
lag = [];

mismatch = 0;


for i = 1:nr
    if mismatch
        % from mismatch to match
        if (btn_visiable(i) == 0 && btn_pressed(i) == 0) || ...
            (btn_visiable(i) == 1 && btn_pressed(i) == 1)
            cur_lag = tsp(i) - start_tsp;
            lag = [lag cur_lag];
            mismatch = 0;
        end
    else
        if btn_visiable(i) ~= btn_pressed(i)
            %(btn_visiable(i) == 1 && btn_pressed(i) == 0) || ...
            %(btn_visiable(i) == 0 && btn_pressed(i) == 1)
            mismatch = 1;
            start_tsp = tsp(i);
        end
    end
end


% statistical features
stats_h = {'sum','mean', 'std', 'q1', 'q3', 'iqr', 'median', 'mode', 'range', ...
    's', 'k', 'mse','En', 'meanTKEO', 'dfa'};
axes_h = {'touch'};
NA = numel(axes_h);
N = numel(stats_h);
A = lag;
for i = 1:NA
    for j = 1:N
        name = strcat(axes_h{i}, '_', stats_h{j});
        header{(i-1)*N + j} = name;
    end
end

Q1 = prctile(A,25);
Q3 = prctile(A,75);
stats_f = [
    sum(A);
    mean(A);
    std(A);
    Q1;
    Q3;
    Q3 - Q1;
    median(A);
    mode(A);
    max(A) - min(A);
    skewness(A);
    kurtosis(A);
    mean(A.*A);
    entropy(A);
    FeatureMeanTKEO(A');
    fastdfa(A');
];

feature = stats_f';

end

