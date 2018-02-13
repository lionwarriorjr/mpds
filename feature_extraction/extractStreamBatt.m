function [ output_args ] = extractStreamBatt( rawFile, featureFile, version, password, fig )
'featureFile', featureFile
% load data in rawFile to workspace
try
    [data] = loadBinAesFile(rawFile, password);
    [r c] = size(data);
    if r > 10
        % extract accel test features
        if strcmp(version, '1')
            [header, feature] = extractStreamBattV1(data);
        end
        % replace ext 'bin' to 'csv'
        csvFile = strrep(featureFile, '.bin.', '.csv.');
        % save features to featureFile
        saveCsvAesFile(header, feature, csvFile, password);
    end
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