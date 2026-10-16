% descriptive script for the methods section
% this script counts the number of:
% 1) discarded epochs 
% 2) interpolated channels
% 3) discarded ICs (adults)
% per participant during artefact rejection

clc
clear
close all

 
%% parameters
par = get_par('experiment', 'newborn');

if strcmp(par.experiment, '2ioi_conga')

    subjects = [2:7 9:15 17:20];

elseif strcmp(par.experiment, 'newborn')

    subjects = [1:3 5:14 16:25 27:29 31:42]; 
end


%% number discarded trials
% retrieve number discarded epochs from the badTrials.csv file
for i = 1:length(subjects)

    fname_bad_trials = fullfile(par.deriv_path, 'preproc', sub_num2str(subjects(i)),...
                 sprintf('%s_badTrials.csv', sub_num2str(subjects(i)))); 

    tbl_bad_trials = readtable(fname_bad_trials); 

    bad_epochs{i} = tbl_bad_trials.trial;

    sum_bad_trial(i) = length(bad_epochs{i});

end

disp(sprintf('experiment: %s', par.experiment))
disp(sprintf('num. discarded trials mean = %.2f std =  %.2f', ...
    mean(sum_bad_trial), std(sum_bad_trial)))
disp(sprintf('range discarded trials = %d -%d', ...
    min(sum_bad_trial), max(sum_bad_trial)))


%% numbers of retained trials

for i_sub = 1:length(subjects)

    sub = subjects(i_sub);
    sub_str = sub_num2str(sub);

    ABprefix = 'AB-100_'; 

    % load data
    if strcmp(par.experiment, 'newborn')

        fpath = fullfile(par.deriv_path, 'preproc', sub_str);
        fname =  sprintf('sub-%03d_task-EEG_filt_ep_%schanint.lw6', ...
            sub, ABprefix);

    else

        fpath = fullfile(par.deriv_path, 'preproc', sub_str);
        fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*.lw6', ...
            sub);

    end

    d = dir(fullfile(fpath, fname)); 

    [header_lf_ep, data_lf_ep] = CLW_load(fullfile(d.folder, d.name));

    sum_kept_trial(i_sub) = length(header_lf_ep.events);

end

disp(sprintf('experiment: %s', par.experiment))
disp(sprintf('num. kept trials mean = %.2f std =  %.2f', ...
    mean(sum_kept_trial), std(sum_kept_trial)))
disp(sprintf('range kept trials = %d -%d', ...
    min(sum_kept_trial), max(sum_kept_trial)))


%%
n_trial_check = (sum_bad_trial+sum_kept_trial)

disp(sprintf('num trial check range range: %d - %d', ...
    min(n_trial_check), max(n_trial_check)))


%% number of interpolated channels
for i = 1:length(subjects)

    fname_bad_chans = fullfile(par.deriv_path, 'preproc', sub_num2str(subjects(i)),...
                 sprintf('%s_badChans.csv', sub_num2str(subjects(i)))); 

    tbl_bad_chans = readtable(fname_bad_chans); 

    bad_chans{i} = tbl_bad_chans.channel;

    sum_bad_chans(i) = length(bad_chans{i});

end


disp(sprintf('experiment: %s', par.experiment))
disp(sprintf('num. interpolated chans mean = %.2f std =  %.2f', ...
    mean(sum_bad_chans), std(sum_bad_chans)))
disp(sprintf('range interpolated channels = %d -%d', ...
    min(sum_bad_chans), max(sum_bad_chans)))


%% retrieve number of channels interpolated by num participants
num_ch_int = unique(sum_bad_chans);

counts = arrayfun(@(x) sum(sum_bad_chans == x), num_ch_int);

[num_ch_int; counts]


%% number discarded ICs 

if strcmp(par.experiment, '2ioi_conga')

    for i = 1:length(subjects)

        fname_bad_ICs = fullfile(par.deriv_path, 'preproc', sub_num2str(subjects(i)),...
            sprintf('%s_badICs.csv', sub_num2str(subjects(i))));

        tbl_bad_ICs = readtable(fname_bad_ICs);

        bad_ICs{i} = tbl_bad_ICs.IC;

        sum_bad_IC(i) = length(bad_ICs{i});

    end

    disp(sprintf('experiment: %s', par.experiment))
    disp(sprintf('num. discarded ICs mean = %.2f std =  %.2f', ...
        mean(sum_bad_chans), std(sum_bad_chans)))
    disp(sprintf('range discarded ICs = %d -%d', ...
        min(sum_bad_chans), max(sum_bad_chans)))

end