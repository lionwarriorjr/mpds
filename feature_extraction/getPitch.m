% @Author: Andong Zhan
% @Date:   2018-02-18 13:40:58
% @Last Modified by:   Andong Zhan
% @Last Modified time: 2018-02-18 14:08:57
% (c) 2018 Max Little. If you use this code, please cite:
% Andong Zhan, et al. "High Frequency Remote Monitoring of Parkinson's Disease
%  via Smartphone: Platform Overview and Medication Response Detection."
% arXiv preprint arXiv:1601.00960 (2016).

function [ pitch ] = getPitch( data, fs, NFFT )
%GETPITCH Summary of this function goes here
%   Calculate pitch
out = pwelch(data, hamming(NFFT),[],NFFT,fs);
[pks locs] = findpeaks(out);
indx = find(pks == max(pks));
indx_max = locs(indx);
pitch = indx_max/length(out)*(fs/2);
end

