% @Author: Andong Zhan
% @Date:   2018-02-18 13:40:58
% @Last Modified by:   Andong Zhan
% @Last Modified time: 2018-02-18 14:06:33
% (c) 2018 Max Little. If you use this code, please cite:
% Andong Zhan, et al. "High Frequency Remote Monitoring of Parkinson's Disease
%  via Smartphone: Platform Overview and Medication Response Detection."
% arXiv preprint arXiv:1601.00960 (2016).

function [header, feature] = extractReactTestV1( data )
% Useage:
% [header, feature] = extractReactTestV1( data )
% Inputs
%    data       - input signal: a matrix where each row is a screen touch
%                 sample. The columns are time, button_visiable, and
%                 button_pressed.
%
% Outputs:
%    header     - the names of the acceleration features
%    feature    - the feature matrix
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

