function [ output_args ] = extractActiveTest(user, session, version, update, password)
%EXTRACT Summary of this function goes here
%   Detailed explanation goes here
% addpath
addpath('InfoTheory');
addpath('mi');
user
session

rootDir = '/home/zad/hopkinspd_server/data/feature';
userDir = fullfile(rootDir, user);

sessionDir = fullfile(userDir, session);
%sessionPatt = strcat('HopkinsPD_', type, '_*', date, '*');

filePatt = strcat('HopkinsPD_test*.temp');

%[userDirs, isdir] = glob(fullfile(rawDir, userPatt));
%nu = size(userDirs);

%for u = 1:nu(1)
%    userDir = userDirs(u);
%    userDir = userDir{1};
%    [sessionDirs, isdir] = glob(fullfile(userDir, sessionPatt));
%    ns = size(sessionDirs);
%    for s = 1:ns(1)
%sessionDir = sessionDirs(s);
%sessionDir = sessionDir{1};
sessionDir
[rawFiles, isdir] = glob(fullfile(sessionDir, filePatt));
aes = 0;
N = size(rawFiles)
for i = 1:N(1)
    rawFile = rawFiles(i);
    rawFile = rawFile{1};
    rawFile


    if ~isempty(strfind(rawFile, 'test1_accel')) || ~isempty(strfind(rawFile, 'test2_accel'))
        featureFile = strrep(rawFile, 'csv.temp', strcat(version, '.csv'));
        if ~update || exist(featureFile, 'file') ~= 2
            extractAccelTest(rawFile, featureFile, version, password,aes);
        end
    elseif ~isempty(strfind(rawFile, 'test0_audio'))
        featureFile = strrep(rawFile, 'raw.temp', strcat(version, '.raw'));
        if ~update || exist(featureFile, 'file') ~= 2
            extractAudioTest(rawFile, featureFile, version, password,aes);
        end
    elseif ~isempty(strfind(rawFile, 'test3_tap'))
        featureFile = strrep(rawFile, '/raw/', '/feature/');
        featureFile = strrep(featureFile, 'csv.temp', strcat(version, '.csv'));
        if ~update || exist(featureFile, 'file') ~= 2
            extractTapTest(rawFile, featureFile, version, password, aes);
        end
    elseif ~isempty(strfind(rawFile, 'test4_react'))
        featureFile = strrep(rawFile, '/raw/', '/feature/');
        featureFile = strrep(featureFile, 'csv.temp', strcat(version, '.csv'));
        if ~update || exist(featureFile, 'file') ~= 2
            extractReactTest(rawFile, featureFile, version, password, aes);
        end
%     elseif strcmp(subtype, 'stream_batt')
%         featureFile = strrep(rawFile, '/raw/', '/feature/');
%         featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
%         if ~update || exist(featureFile, 'file') ~= 2
%             extractStreamBatt(rawFile, featureFile, version, password);
%         end
%     elseif strcmp(subtype, 'stream_meta')
%         featureFile = strrep(rawFile, '/raw/', '/feature/');
%         featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
%         if ~update || exist(featureFile, 'file') ~= 2
%             extractStreamMeta(rawFile, featureFile, version, password);
%         end
%     elseif strcmp(subtype, 'stream_calllog')
%         featureFile = strrep(rawFile, '/raw/', '/feature/');
%         featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
%         if ~update || exist(featureFile, 'file') ~= 2
%             extractStreamCalllog(rawFile, featureFile, version, password);
%         end
%     elseif strcmp(subtype, 'stream_smslog')
%         featureFile = strrep(rawFile, '/raw/', '/feature/');
%         featureFile = strrep(featureFile, 'csv', strcat(version, '.csv'));
%         if ~update || exist(featureFile, 'file') ~= 2
%             extractStreamSmslog(rawFile, featureFile, version, password);
%         end
    end
end
%end
%end
%end

