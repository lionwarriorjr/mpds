function [ output_args ] = extractTapTest( rawFile, featureFile, version, password, aes )
%EXTRACTACCELTEST Summary of this function goes here
%   Detailed explanation goes here
% 'rawFile', rawFile
'featureFile', featureFile
% 'Version', version

% load data in rawFile to workspace
readHeader = 0;
try
    data = loadAesFile(rawFile, password, readHeader, aes);

    %% preprocessing
    [nr, nc] = size(data);

    % make sure the time vector is increasing

    if ~all(diff(data(:,1))>0)
        data = sortrows(data, 1);
        % delete redundant sample
        for i = 1:nr-1
            j = nr-i;
            if data(j,1) == data(j+1,1)
                data(j+1,:) = [];
            end
        end
    end

    [nr, nc] = size(data);

    % extract accel test features
    if strcmp(version, '1')
        [header, feature] = extractTapTestV1(data);
    elseif strcmp(version, '2')
        [header, feature] = extractTapTestV2(data);
    end

    % save features to featureFile
    saveCsvAesFile(header, feature, featureFile, password);
catch ME
    eid = ME.identifier
    if ~(strcmp(eid, 'MATLAB:csvread:FileNotFound') || ...
            strcmp(eid, 'MATLAB:textscan:EmptyFormatString')...
            )
        featureFile
        ME
    end
end
end

