function [ output_args ] = extractStream(userPatt, date, type, subtype, version, update, fig, password )
%EXTRACT Summary of this function goes here
%   Detailed explanation goes here
% addpath
addpath('InfoTheory');
addpath('mi');

rootDir = '/home/zad/hopkinspd_server/data';
inDir = fullfile(rootDir, 'incoming/');

sessionPatt = strcat('HopkinsPD_', type, '_*', date, '*.zip.aes');

filePatt = strcat('HopkinsPD_', subtype, '_*');

[userDirs, isdir] = glob(fullfile(inDir, userPatt));
nu = size(userDirs);

for u = 1:nu(1)
    userDir = userDirs(u);
    userDir = userDir{1};
    [sessionDirs, isdir] = glob(fullfile(userDir, sessionPatt));
    ns = size(sessionDirs);
    for s = 1:ns(1)
        sessionDir = sessionDirs(s);
        sessionDir = sessionDir{1};
        try
            % decrypt and unzip it
            descryptAndUnzip(sessionDir,password);
            %%%%%%
            [rawFiles, isdir] = glob(fullfile('temp', filePatt));
            N = size(rawFiles);
            for i = 1:N(1)
                rawFile = rawFiles(i);
                rawFile = rawFile{1};
                if strcmp(subtype, 'stream_accel')
                    featureFile = strrep(rawFile, 'temp', sessionDir(1:length(sessionDir)-8));
                    featureFile = strrep(featureFile, '.raw', strcat(version, '.csv'));
                    featureFile = strrep(featureFile, '/incoming/', '/feature/');
                    if ~update || exist(featureFile, 'file') ~= 2
                        extractStreamAccel(rawFile, featureFile, version, password, fig);
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
        catch ME
            ME
        end
    end
end
end

