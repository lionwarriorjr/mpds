function [ output_args ] = extractAudioTest( rawFile, featureFile, version, password, aes )
%EXTRACTACCELTEST Summary of this function goes here
%   Detailed explanation goes here
% 'rawFile', rawFile
'featureFile', featureFile
% 'Version', version

% load data in rawFile to workspace
try
    [data, fs] = loadAudioAesFile(rawFile, password, aes);


    % extract accel test features
    if strcmp(version, '1')
        [header, feature] = extractAudioTestV1(data, fs);
    end
    % replace ext 'raw' to 'csv'
    csvFile = strrep(featureFile, '.raw', '.csv');
    % save features to featureFile
    saveCsvAesFile(header, feature, csvFile, password);
catch ME
    eid = ME.identifier;
    if ~(strcmp(eid, 'MATLAB:csvread:FileNotFound') || ...
            strcmp(eid, 'MATLAB:textscan:EmptyFormatString')...
            )
        featureFile
        ME
    end
end
end

