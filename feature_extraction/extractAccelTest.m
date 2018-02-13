function [ output_args ] = extractAccelTest( rawFile, featureFile, version, password, aes )
%EXTRACTACCELTEST Summary of this function goes here
%   Detailed explanation goes here
'rawFile', rawFile
'featureFile', featureFile
% 'Version', version

% load data in rawFile to workspace
try
    data = loadAesFile(rawFile, password, 0 , aes);

    %% preprocessing
    [nr, nc] = size(data);
    tsp = data(1,1);

    % make sure the time vector is increasing
    invalid_idx = 0;

    for i = 2:nr
        if data(i,1) <= tsp
            invalid_idx = 1;
            break;
        else
            tsp = data(i,1);
        end
    end

    if invalid_idx == 1
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

    % filter first 25% and last 25% data
    q = round(nr/4);
    data = data(q:nr-q, :);

    %% extract accel test features
    if strcmp(version, '1')
        [header, feature] = extractAccelTestV1(data);
    elseif strcmp(version, '2')
        [header, feature] = extractAccelTestV2(data);
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

