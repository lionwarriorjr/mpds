function [header, feature] = extractTapTestV2( data )
%EXTRACTACCELTESTV2 Summary of this function goes here
%   Detailed explanation goes here
[nr, nc] = size(data);
tsp = data(:,1) - data(1,1);
x = data(:,2);
y = data(:,3);

maxX = max(x);
minX = min(x);
mid = (maxX - minX)/2 + minX;
s = 0;
start = 0;
prev = 0;
j = 1;
for i = 1:nr
    if s == 0
        % init
        s = sign(x(i) - mid);
        start = tsp(i);
    else
        new_sign = sign(x(i) - mid);
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

ntaps = length(Tstay);
% scale coordinates from [0, roundn(max(x)] to [0,1]
scaler = roundn(max(x),2);
scaleX = x/scaler;
scaleY = y/scaler;

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
% cross features on stay and move axis
cross_h = {'t_xcorr','t_mi','t_xEn'};

header = [header cross_h];

r = corrcoef(A(:,1),A(:,2));
corr = r(1,2);
MI = mi(A(:,1),A(:,2));
xEn = entropy(A(:,1)) + ...
    relativeEntropy(round(A(:,1).*10000),round(A(:,2).*10000));

cross_f = [corr;MI;xEn];


feature = [stats_f cross_f'];


A = [scaleX, scaleY];
stats_h = {'mean', 'std', 'q1', 'q3', 'iqr', 'median', 'mode', 'range', ...
    's', 'k', 'mse','En', 'zcr', 'dfc', 'dfc_amp', 'meanTKEO', 'ar1', 'dfa'};
axes_h = {'scaled_x','scaled_y'};
NA = numel(axes_h);
N = numel(stats_h);
for i = 1:NA
    for j = 1:N
        name = strcat(axes_h{i}, '_', stats_h{j});
        header2{(i-1)*N + j} = name;
    end
end
% dominant frequency component
%F = 1:10;
minF = 0.5;
maxF = 20;
[pxx, f] = plomb(A,tsp,maxF);
ind = f > minF;
f = f(ind);
pxx = pxx(ind, :);


[amp, maxi] = max(pxx);
dfc = f(maxi)';
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
    entropy(A(:,1)) entropy(A(:,2)) ;
    ZCR(normalize(A(:,1))) ZCR(normalize(A(:,2))) ;
    dfc;
    amp;
    FeatureMeanTKEO(A(:,1)) FeatureMeanTKEO(A(:,2));
    FeatureAR1(A(:,1)) FeatureAR1(A(:,2));
    fastdfa(A(:,1)) fastdfa(A(:,2));
];
stats_f = reshape(stats_f, [1,N*2]);
% cross features on x,y axis
cross_h = {'scaled_xy_xcorr','scaled_xy_mi','scaled_xy_xEn'};

header = [{'ntaps'} header header2 cross_h];

r = corrcoef(A(:,1),A(:,2));
corr = r(1,2);
MI = mi(A(:,1),A(:,2));
xEn = entropy(A(:,1)) + ...
    relativeEntropy(round(A(:,1).*10000),round(A(:,2).*10000));

cross_f = [corr;MI;xEn];


feature = [[ntaps] feature stats_f cross_f'];

% add features about left and right touch locations
left_idx = x < mid;
right_idx = x > mid;
scaleL = [scaleX(left_idx) scaleY(left_idx)];
tspL = tsp(left_idx);
scaleR = [scaleX(right_idx) scaleY(right_idx)];
tspR = tsp(right_idx);
[scaleL_h, scaleL_f] = getStatFeatures(scaleL, tspL, {'scaled_x_left', 'scaled_y_left'});
[scaleR_h, scaleR_f] = getStatFeatures(scaleR, tspR, {'scaled_x_right', 'scaled_y_right'});
header = [header scaleL_h scaleR_h];
feature = [feature scaleL_f, scaleR_f];

% add cross fetures on scaled_xy_left and scaled_xy_right
[cross_h, cross_f] = getCrossFeatures(scaleL, 'scaled_xy_left');
header = [header cross_h];
feature = [feature cross_f];
[cross_h, cross_f] = getCrossFeatures(scaleR, 'scaled_xy_right');
header = [header cross_h];
feature = [feature cross_f];

end

