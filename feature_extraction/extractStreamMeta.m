function [ output_args ] = extractStreamMeta( rawFile, featureFile, version, password, fig )
'featureFile', featureFile
% load data in rawFile to workspace
try
    cmd = sprintf('ln -s %s %s', rawFile, featureFile);
    system(cmd);

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