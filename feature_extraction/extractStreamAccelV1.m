function [header, feature] = extractStreamAccelV1( data, fig )


[nr, nc] = size(data);
header = {};
X = zeros(nr, nc+4);
X(:,1:4) = data;
[a,e,r] = cart2sph(data(:,2),data(:,3),data(:,4));
X(:,5) = a;
X(:,6) = e;
X(:,7) = r;
T = X(:,1);
A = X(:,2:7);
% statistical features
stats_h = {'mean', 'std', 'q1', 'q3', 'iqr', 'median', 'mode', 'range', ...
    's', 'k', 'mse','En', 'zcr', 'dfc', 'dfc_amp', 'meanTKEO', 'ar1', 'dfa'};
axes_h = {'x','y','z','a','e','r'};
NA = numel(axes_h);
N = numel(stats_h);
for i = 1:NA
    for j = 1:N
        name = strcat(axes_h{i}, '_', stats_h{j});
        header{(i-1)*N + j} = name;
    end
end
% dominant frequency component
%F = 1:10;

minF = 0;
maxF = 20;

% use plomb but too slow
%[pxx, f] = plomb(A,T,maxF);
%ind = f > minF;
%f = f(ind);
%pxx = pxx(ind, :);
%[amp, maxi] = max(pxx);
%dfc = f(maxi)';
Fs = nr/data(nr,1);
dfc = zeros(1,6);
amp = zeros(1,6);
for i = 1:6
    [d, a] = DFCfft(A(:,i)', Fs, minF, maxF);
    dfc(i) = d;
    amp(i) = a;
end

Q1 = prctile(A,25);
Q3 = prctile(A,75);
meanA = A;
m = mean(A);
for i = 1:nr
    meanA(i,:) = meanA(i,:) - m;
end
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
    entropy(A(:,1)) entropy(A(:,2)) entropy(A(:,3)) entropy(A(:,4))...
    entropy(A(:,5)) entropy(A(:,6));
    ZCR(meanA(:,1)) ZCR(meanA(:,2)) ZCR(meanA(:,3)) ZCR(meanA(:,4))...
    ZCR(meanA(:,5)) ZCR(meanA(:,6));
    dfc;
    amp;
    FeatureMeanTKEO(A(:,1)) FeatureMeanTKEO(A(:,2))...
    FeatureMeanTKEO(A(:,3)) FeatureMeanTKEO(A(:,4))...
    FeatureMeanTKEO(A(:,5)) FeatureMeanTKEO(A(:,6));
    FeatureAR1(A(:,1)) FeatureAR1(A(:,2)) FeatureAR1(A(:,3)) FeatureAR1(A(:,4)) ...
    FeatureAR1(A(:,5)) FeatureAR1(A(:,6));
    fastdfa(A(:,1)) fastdfa(A(:,2)) fastdfa(A(:,3)) fastdfa(A(:,4))...
    fastdfa(A(:,5)) fastdfa(A(:,6));
];
stats_f = reshape(stats_f, [1,N*6]);
% cross features on x,y,z axis
cross_h = {'xcorr','mi','xEn'};
corr = zeros(1,6);
MI = zeros(1,6);
xEn = zeros(1,6);
idx = 1;
nh = numel(header);
for i = 1:2
    for j = i+1:3
        for k = 1:numel(cross_h)
            header{nh + (idx-1)*numel(cross_h) + k} ...
                = strcat(axes_h{i}, '_',axes_h{j}, '_', cross_h{k});
        end
        r = corrcoef(A(:,i),A(:,j));
        corr(1,idx) = r(1,2);
        MI(1,idx) = mi(A(:,i),A(:,j));
        xEn(1,idx) = entropy(A(:,i)) + ...
            relativeEntropy(round(A(:,i).*10000),round(A(:,j).*10000));
        idx = idx + 1;
    end
end
% cross features on a,e,r
nh = numel(header);
nidx = idx;
for i = 4:5
    for j = i+1:6
        for k = 1:numel(cross_h)
            header{nh + (idx-nidx)*numel(cross_h) + k} ...
                = strcat(axes_h{i}, '_',axes_h{j}, '_', cross_h{k});
        end
        r = corrcoef(A(:,i),A(:,j));
        corr(1,idx) = r(1,2);
        MI(1,idx) = mi(A(:,i),A(:,j));
        xEn(1,idx) = entropy(A(:,i)) + ...
            relativeEntropy(round(A(:,i).*10000),round(A(:,j).*10000));
        idx = idx + 1;
    end
end


cross_f = [corr;MI;xEn];
cross_f = reshape(cross_f, 1, 6*numel(cross_h));
feature = [stats_f cross_f];


end

