function [header, feature] = extractTapTestV1( data )
%EXTRACTACCELTESTV1 Summary of this function goes here
%   Detailed explanation goes here
[nr, nc] = size(data);
tsp = data(:,1);
y = data(:,2);
maxY = max(y);
minY = min(y);
mid = (maxY - minY)/2 + minY;
s = 0;
start = 0;
prev = 0;
j = 1;
for i = 1:nr
    if s == 0
        % init
        s = sign(y(i) - mid);
        start = tsp(i);
    else
        new_sign = sign(y(i) - mid);
        if new_sign ~= s
            Tstay(j) = prev - start;
            Tmove(j) = tsp(i) - prev;
            j = j+1;
            start = tsp(i);
            s = new_sign;
        end
    end
    prev = tsp(i);
end



% statistical features
stats_h = {'mean', 'std', 'q1', 'q3', 'iqr', 'median', 'mode', 'range', ...
    's', 'k', 'mse','En', 'meanTKEO', 'ar1', 'dfa'};
axes_h = {'stay','move'};
NA = numel(axes_h);
N = numel(stats_h);
A = [Tstay;Tmove]';
for i = 1:NA
    for j = 1:N
        name = strcat(axes_h{i}, '_', stats_h{j});
        header{(i-1)*N + j} = name;
    end
end

Q1 = prctile(A,25);
Q3 = prctile(A,75);
stats_f = [
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
    entropy(A(:,1)) entropy(A(:,2));
    FeatureMeanTKEO(A(:,1)) FeatureMeanTKEO(A(:,2));
    FeatureAR1(A(:,1)) FeatureAR1(A(:,2));
    fastdfa(A(:,1)) fastdfa(A(:,2));
];
stats_f = reshape(stats_f, [1,N*2]);
% cross features on x,y,z and acc axis
cross_h = {'xcorr','mi','xEn'};

header = [header cross_h];

r = corrcoef(A(:,1),A(:,2));
corr = r(1,2);
MI = mi(A(:,1),A(:,2));
xEn = entropy(A(:,1)) + ...
    relativeEntropy(round(A(:,1).*10000),round(A(:,2).*10000));

cross_f = [corr;MI;xEn];
feature = [stats_f cross_f'];


end

