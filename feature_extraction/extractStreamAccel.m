function [ output_args ] = extractStreamAccel( rawFile, featureFile, version, password, fig )
%EXTRACTACCELTEST Summary of this function goes here
%   Detailed explanation goes here
'rawFile', rawFile
'featureFile', featureFile
% 'Version', version

% load data in rawFile to workspace
H = [];
F = [];
try
    data = loadRawFile(rawFile);

    %% preprocessing

    % change time diff to tsp
    [nr, nc] = size(data);
    pre = data(1,1);
    for i = 2:nr
        data(i,1) = pre + data(i,1);
        pre = data(i,1);
    end
    % sort rows based on the tsp
    data = sortrows(data, 1);

    % delete redundant sample (only for plomb)
    %for i = 1:nr-1
    %    j = nr-i;
    %    if data(j,1) == data(j+1,1)
    %        data(j+1,:) = [];
    %    end
    %end
    [nr, nc] = size(data);

    % average sampling hz
    hz = nr/data(nr,1);
    min_hz = 20;

    % slide data into windows with window size and slide in seconds
    win = 2*60;
    slide = 1*60;

    % todo while loop to find each window
    win_start = data(1,1);
    win_end = win_start+win;
    while win_start < data(nr,1)
        win_start
        win_data = data(data(:,1)>=win_start & data(:,1)<win_end,:);
        % todo check there are enough samples in each window
        [ns ig] = size(win_data);
        if ns > min_hz*win
            % make tsp starts from zero
            win_data(:,1) = win_data(:,1) - win_data(1,1);
            % todo generate one feature vector for each window
            %% extract accel test features
            if strcmp(version, '1')
                [header, feature] = extractStreamAccelV1(win_data, fig);
            end
            if isempty(H)
                H = ['sec' header];
            end
            if isempty(F)
                F = [win_start feature];
            else
                F = [F;win_start feature];
            end
        end
        win_start = win_start + slide;
        win_end = win_end + slide;
    end
    if isempty(H)
        disp('empty data')
    else
        % save features to featureFile
        saveCsvAesFile(H, F, featureFile, password);
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

