clc; clear; close all;

channels = 22;
fs = 250;

t_min = 2.0;
t_max = 6.0;

type = 2;

raw_data_path = 'E:\Polyspace\BCI_Proj\BCICIV2a\';
%{
file_and_folders = dir(raw_data_path);
file_names = {file_and_folders(~[file_and_folders.isdir]).name};

for i = 1 : numel(file_names)
    file_path = strcat(raw_data_path, file_names{i});
    disp(file_path)
end
%}
save_path = 'E:\Polyspace\BCI_Proj\BCICIV2a\lowpass\';

for i = 1 : 9
    T_file_path = strcat(raw_data_path, 'A0', num2str(i), 'T.mat');
    E_file_path = strcat(raw_data_path, 'A0', num2str(i), 'E.mat');
    
    subj_save_path = strcat(save_path, 'A0', num2str(i), '.mat');

    [T_trials, T_labels] = filter_and_crop(T_file_path, ...
        t_min, t_max, type, fs, channels);
    [E_trials, E_labels] = filter_and_crop(E_file_path, ...
        t_min, t_max, type, fs, channels);

    save(subj_save_path, 'T_trials', 'T_labels', 'E_trials', 'E_labels')
end

save(strcat(save_path, 'info.mat'), 'channels', 'fs', 't_min', 't_max', 'type')

function [trials, labels] = filter_and_crop(path, t_min, t_max, type, fs, channels)
    % a version for band pass filtering
    fc = 250; % sampling rate
    Wl = 4; Wh = 40; % pass band
    Wn = [Wl*2 Wh*2]/fc;
    [b,a]=cheby2(6,60,Wn);

    % Load data
    data = load(path).data;
    n_runs = numel(data);

    n_valid_runs = 0; % to record the number of valid runs
    trials = zeros(48 * n_runs, channels, fs * (t_max - t_min));
    labels = zeros(48 * n_runs, 1);

    for j = 1 : n_runs
        % read data
        run = data{j};
        X = run.X(:, 1:channels);
        events = run.trial;
        y = run.y;
        
        % dismiss idle runs
        if size(events, 1) == 0
            continue
        else
            n_valid_runs = n_valid_runs + 1;
        end

%{
        % filter EEG signals with band [Wl, Wh]
        if type == 0
            filtered_X = X;
        elseif type == 1
            filtered_X = bandpass(X, [4.0, 40.0], 250.0);
        elseif type == 2
            filtered_X = lowpass(X, 40.0, 250.0);
        else
            error('Select filtering type from none (0), bandpass (1) or lowpass (2)')
        end
%}

        % Crop signals between [t_min, t_max] after each event
        for k = 1 : 48
            start_point = events(k) + fs * t_min;
            end_point = events(k) + fs * t_max - 1;
            cropped_trial = X(start_point : end_point, :);
            % disp(size(cropped_trial))
            % filter EEG signals with band [Wl, Wh]
            if type == 0
                filtered_trial = cropped_trial;
            elseif type == 1
                filtered_trial = filtfilt(b, a, cropped_trial);
            elseif type == 2
                filtered_trial = lowpass(cropped_trial, 40.0, 250.0);
            else
                error('Select filtering type from none (0), bandpass (1) or lowpass (2)')
            end
            trials(48*(n_valid_runs-1)+k, :, :) = filtered_trial';
        end

        % record labels
        labels(48*(n_valid_runs-1)+1 : 48*n_valid_runs) = y;
    end

    trials = trials(1:288, :, :);
    labels = labels(1:288);
    
end