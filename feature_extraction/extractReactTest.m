function [ output_args ] = extractTapTest( rawFile, featureFile, version, password, aes )
%EXTRACTACCELTEST Summary of this function goes here
%   Detailed explanation goes here
% 'rawFile', rawFile
'featureFile', featureFile
% 'Version', version

% load data in rawFile to workspace
try
    data = loadAesFile(rawFile, password,0, aes);
    % extract accel test features
    if strcmp(version, '1')
        [header, feature] = extractReactTestV1(data);
    end

    % save features to featureFile
    saveCsvAesFile(header, feature, featureFile, password);
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

