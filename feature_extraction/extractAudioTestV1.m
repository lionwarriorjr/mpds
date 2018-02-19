% @Author: Andong Zhan
% @Date:   2018-02-18 13:40:58
% @Last Modified by:   Andong Zhan
% @Last Modified time: 2018-02-18 14:04:10
% (c) 2018 Max Little. If you use this code, please cite:
% Andong Zhan, et al. "High Frequency Remote Monitoring of Parkinson's Disease
%  via Smartphone: Platform Overview and Medication Response Detection."
% arXiv preprint arXiv:1601.00960 (2016).

function [header, feature] = extractAudioTestV1( data, fs )
% Useage:
% [header, feature] = extractAudioTestV1( data )
% Inputs
%    data       - input signal: audio signal.
%    fs         - audio frequency: e.g., 44100
% Outputs:
%    header     - the names of the acceleration features
%    feature    - the feature matrix
%


% set frame size to 0.5 second
frameSizeInSecs = 0.5;
frameSize = fs*frameSizeInSecs;
nframes = size(data)/frameSize;
nframes = floor(nframes(1));
%% find the seconds have voice in the sample, i.e., volume > 0.1
frames = cell(1,nframes);
amp = zeros(1,nframes);
pitch = zeros(1,nframes);
num = 1;
startF = 0;
for i = 1:nframes
    frames{i} = data((1+(i-1)*frameSize):i*frameSize,1);
    amp(i) = max(abs(frames{i}));
end

% amplitude threshold
ampT = (quantile(amp,0.9) - min(amp))/4 + min(amp);
for i = 1:nframes
    if amp(i) > ampT
        if startF == 0
            startF = i;
            voice(num).start = i;
        end
    else
        if startF ~=0
            voice(num).end = i;
            voice(num).len = (i - voice(num).start) * frameSizeInSecs;
            num = num + 1;
            startF = 0;
        end
    end
end

if startF ~= 0
    voice(num).end = nframes;
    voice(num).len = (nframes - voice(num).start + 1) * frameSizeInSecs;
end




%% select the longest voice V to process
len = 0;
s = size(voice);
for i = 1:s(2)
    if len < voice(i).len
        len = voice(i).len;
        V = voice(i);
    end
end
% pick the middle part of voice
range = V.end - V.start + 1;
V.start = V.start + floor(range/4);
V.end = V.end - floor(range/4);

%% voice feature extraction
% V.len: voice length in seconds
% voice mean amplitude
fprintf('voice length: %d', V.len);
V.amp_mean = mean(amp(V.start:V.end));
V.amp_std = std(amp(V.start:V.end));
V.data = data((1+(V.start-1)*frameSize):V.end*frameSize,1);
V.amp_p1 = polyfit(V.start:V.end, amp(V.start:V.end), 1);
V.amp_p2 = polyfit(V.start:V.end, amp(V.start:V.end), 2);
V.amp_dfa = fastdfa(amp(V.start:V.end)');
% V.pitch: dorminant voice frequency in total

for i = V.start:V.end
    pitch(i) = getPitch(frames{i}, fs, 10240);
end

% TODO pitch, amplitude, variance & trend

V.pitch = getPitch(V.data, fs, 10240);
V.pitch_std = std(pitch(V.start:V.end));
V.pitch_p1 = polyfit(V.start:V.end, pitch(V.start:V.end), 1);
V.pitch_p2 = polyfit(V.start:V.end, pitch(V.start:V.end), 2);
V.pitch_dfa = fastdfa(pitch(V.start:V.end)');


% statistical features
header = {'length','amp_mean','amp_std', 'amp_p1_1', 'amp_p1_0', ...
    'amp_p2_2','amp_p2_1','amp_p2_0', 'amp_dfa',...
    'pitch','pitch_std','pitch_p1_1','pitch_p1_0','pitch_p2_2',...
    'pitch_p2_1','pitch_p2_0','pitch_dfa'};

feature = [V.len V.amp_mean V.amp_std V.amp_p1 V.amp_p2 V.amp_dfa...
    V.pitch V.pitch_std V.pitch_p1 V.pitch_p2 V.pitch_dfa];




end

