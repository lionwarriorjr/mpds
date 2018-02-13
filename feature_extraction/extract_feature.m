function [ output_args ] = extract_feature(userPatt, type, subtype, version, update, fig )
%EXTRACT Summary of this function goes here
%   Detailed explanation goes here
% addpath
addpath('InfoTheory');
addpath('mi');

rootDir = '/home/zad/hopkinspd_server/data';
rawDir = fullfile(rootDir, 'raw/');

sessionPatt = strcat('HopkinsPD_', type, '_*');

filePatt = strcat('HopkinsPD_', type, '*', subtype, '*');

[userDirs, isdir] = glob(fullfile(rawDir, userPatt));
password = -1;
nu = size(userDirs);
aes = 0; % raw file has already been decrypted
for u = 1:nu(1)
    userDir = userDirs(u);
    userDir = userDir{1}
    [sessionDirs, isdir] = glob(fullfile(userDir, sessionPatt));
    ns = size(sessionDirs);
    for s = 1:ns(1)
        sessionDir = sessionDirs(s);
        sessionDir = sessionDir{1}
        [rawFiles, isdir] = glob(fullfile(sessionDir, filePatt));
        N = size(rawFiles);
        for i = 1:N(1)
            rawFile = rawFiles(i);
            rawFile = rawFile{1}
            if strcmp(subtype, 'accel')
                featureFile = strrep(rawFile, '/raw/', '/feature/');
                featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
                if ~update || exist(featureFile, 'file') ~= 2
                    extractAccelTest(rawFile, featureFile, version, password, aes);
                end
            elseif strcmp(subtype, 'audio')
                featureFile = strrep(rawFile, '/raw/', '/feature/');
                featureFile = strrep(featureFile, 'raw', strcat(version, '.csv'));
                if ~update || exist(featureFile, 'file') ~= 2
                    extractAudioTest(rawFile, featureFile, version, password, aes);
                end
            elseif strcmp(subtype, 'tap')
                featureFile = strrep(rawFile, '/raw/', '/feature/');
                featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
                if ~update || exist(featureFile, 'file') ~= 2
                    extractTapTest(rawFile, featureFile, version, password, aes);
                end
            elseif strcmp(subtype, 'react')
                featureFile = strrep(rawFile, '/raw/', '/feature/');
                featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
                if ~update || exist(featureFile, 'file') ~= 2
                    extractReactTest(rawFile, featureFile, version, password, aes);
                end
            elseif strcmp(subtype, 'stream_batt')
                featureFile = strrep(rawFile, '/raw/', '/feature/');
                featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
                if ~update || exist(featureFile, 'file') ~= 2
                    extractStreamBatt(rawFile, featureFile, version, password);
                end
            elseif strcmp(subtype, 'stream_meta')
                featureFile = strrep(rawFile, '/raw/', '/feature/');
                featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
                if ~update || exist(featureFile, 'file') ~= 2
                    extractStreamMeta(rawFile, featureFile, version, password);
                end
            elseif strcmp(subtype, 'stream_calllog')
                featureFile = strrep(rawFile, '/raw/', '/feature/');
                featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
                if ~update || exist(featureFile, 'file') ~= 2
                    extractStreamCalllog(rawFile, featureFile, version, password);
                end
            elseif strcmp(subtype, 'stream_smslog')
                featureFile = strrep(rawFile, '/raw/', '/feature/');
                featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
                if ~update || exist(featureFile, 'file') ~= 2
                    extractStreamSmslog(rawFile, featureFile, version, password);
                end
            end
        end
    end
end
end

