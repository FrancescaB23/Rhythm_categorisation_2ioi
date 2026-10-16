% Script to compute grand averaging selecting epochs based on zSNR

clc
clear
close all


%% parameters
par = get_par('experiment', 'newborn'); 

subjects = [num2cell([1:3 5:14 16:25 27:29 31:42])];

harmonics = [1:6];

roi_name = 'front'; 

z_threshold = norminv(1 - .05); 
subj_grand_avg_suffix = sprintf('grand-%d-zSNR',length(subjects));


%% select trials from individual participants based on zSNR
% frequencies of interest
frex = 1 / par.cycle_dur * harmonics;

% condition codes
labels = arrayfun(@(k) sprintf('S%03d', k), 1:par.n_cond, ...
    'UniformOutput', false);

% initialise structures 
data_lf_all =[];
events_lf_all = [];
header_template = [];
count_trials  = zeros(par.n_cond, length(subjects));
subj_no_epochs = {};


for i_sub=1:length(subjects)

    sub = subjects{i_sub};

    sub_str = sub_num2str(sub)


    % load data
    fpath = fullfile(par.deriv_path, 'preproc', sub_str);
    fname =  sprintf('sub-%03d_task-EEG_filt_ep_AB-100_chanint.lw6', sub);
    
    d = dir(fullfile(fpath, fname));

    [header_lf_ep, data_lf_ep] = CLW_load(fullfile(d.folder, d.name));

    % remove some more 'external' channels for consistency with Edalati
    % et al. 2023
    if strcmp(par.experiment, 'newborn')

        chans_outer_ring = {'E14' 'E17' 'E21' 'E43' 'E48' 'E38' 'E44' ...
            'E49' 'E56' 'E107' 'E113' 'E120' 'E1' 'E32' 'E25' 'E8' ...
            'E121' 'E114' 'E99' 'E94' 'E88' 'E81' 'E73' 'E68' 'E63' ...
            'E119' 'E125' 'E126' 'E127' 'E128'};

        chans_outer_ring_idx = find(ismember( ...
            {header_lf_ep.chanlocs.labels}, chans_outer_ring));

        data_lf_ep(:,chans_outer_ring_idx,:,:,:,:) = [];

        header_lf_ep.datasize = size(data_lf_ep);
        header_lf_ep.chanlocs(chans_outer_ring_idx) = [];
    end

    % re-segment
    codes = unique({header_lf_ep.events.code});

    [header_lf_ep, data_lf_ep] = segment_safe(...
        header_lf_ep, ...
        data_lf_ep, ...
        codes,...
        'x_start', 0, ...
        'x_duration', par.trial_dur);

    % downsample
    [header_lf_ep, data_lf_ep] = RLW_downsample(header_lf_ep, data_lf_ep, ...
        'x_downsample_ratio', par.decim_factor);

    % rereference
    [header_lf_ep, data_lf_ep] = RLW_rereference(header_lf_ep, data_lf_ep, ...
        'apply_list', {header_lf_ep.chanlocs.labels}, ...
        'reference_list', par.eeg_ref_chans);
        


    % select and average channels of interest
    if strcmp(roi_name, 'front')

        roi_chans = {
            'E19'
            'E11'
            'E4'
            'E20'
            'E12'
            'E5'
            'E118'
            'E30'
            'E13'
            'E6'
            'E112'
            'E105'
            'E7'
            'E106'
            'E29'
            'E111'
            'E24'
            'E124'
            'E28'
            'E117'
            };

    end

    [header_ch, data_ch] = RLW_arrange_channels( ...
        header_lf_ep, data_lf_ep, roi_chans);

    [header_ch, data_ch] = RLW_pool_channels( ...
        header_ch, data_ch, roi_chans, ...
        'mixed_channel_label', roi_name, ...
        'keep_original_channels', false);


    % compute amplitude FFT
    [header_mX, data_mX] = RLW_FFT(header_ch, data_ch);

    fmax = (1/header_ch.xstep) / 2;
    [header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax);
    data_mX(:, :, :, :, :, 1) = 0;
    header_mX.name = sprintf('%s_snr-none_mX', sub_str);


    % frequency vector
    freq = [0 : header_mX.datasize(end)-1] * header_mX.xstep;

    % zSNR
    z_snr_roi = get_z_snr(data_mX, freq, frex, par.snr_bins(1), par.snr_bins(2));

    % select subset of trials depending on their zSNR > threshold
    idx_zsnr = z_snr_roi > z_threshold;

    if any(idx_zsnr)

        count_trials(:, i_sub) = countcats(categorical( ...
            {header_ch.events(idx_zsnr).code}, labels));

        % select epochs based on zSNR
        % {header_ch.events(idx_zsnr).code}
        [header_sel, data_sel] = RLW_arrange_epochs( ...
            header_lf_ep, data_lf_ep, find(idx_zsnr));

        % concatenate individual data and header
        data_lf_all = cat(1, data_lf_all, data_sel);
        events_lf_all = [events_lf_all, header_sel.events];

        header_template = header_sel;


    else
        warning('%s: no epochs with zSNR > %.2f, subject skipped', ...
            sub_str, z_threshold);
        subj_no_epochs{end+1} = sub_str;
    end


end

fprintf('%d/%d subjects contributed epochs.\n', ...
    length(subjects) - numel(subj_no_epochs), length(subjects));
if ~isempty(subj_no_epochs)
    fprintf('Subjects without selected epochs: %s\n', strjoin(subj_no_epochs, ', '));
end


%% prepare grand average data
% adjust header information to reflect the data
ep_fix = num2cell(1:numel(events_lf_all));
[events_lf_all.epoch] = ep_fix{:};

header_lf_all = header_template;

header_lf_all.name = ['sub-', subj_grand_avg_suffix, header_template.name(8:end)];
header_lf_all.datasize = size(data_lf_all);
header_lf_all.events = events_lf_all;

num_trials_by_cond = sum(count_trials, 2)

% warn if a condition ended up with no trials at all
empty_cond = find(num_trials_by_cond == 0);
if ~isempty(empty_cond)
    warning('No trials selected for condition(s): %s', ...
        strjoin(labels(empty_cond), ', '));
end

% average trials by condition
[header_lf_avg, data_lf_avg] = average_trials_per_cond( ...
    header_lf_all, data_lf_all, par, 'concatenate_epochs', false);

header_lf_avg.name = sprintf('%s_rr-m', header_lf_avg.name);


%% save grand average file

save_path = fullfile(par.deriv_path, 'preproc', ...
    sprintf('sub-%s', subj_grand_avg_suffix));

if ~exist(save_path)
    mkdir(save_path);
end

CLW_save(save_path, header_lf_avg, data_lf_avg);

% save participants list
participant_list_filename = sprintf('sub-%s_participants_list.txt', ...
    subj_grand_avg_suffix);
fid = fopen(fullfile(save_path, participant_list_filename), 'w');  % Open file for writing
fprintf(fid, 'participants\n');      % Write header
fprintf(fid, '%d\n', subjects{:});             % Write each element on a new line
fclose(fid);

