% compute ERP of grand average selected newborn data
% It includes an extra step of artefact detection based on a threshold and
% it performs equalisation of number of epochs across conditions

clc
close all
clear


%% parameters
par = get_par('experiment', 'newborn');
subjects = [1 5 7 8 9 10 12 17 18 22 23 27 28 31 33 35 36 37 38];

roi_name = 'front';

threshold_range_chunk = 100; % threshold for range for chunk


%%
for i_sub = 1:length(subjects)


    sub = subjects(i_sub)


    %% prepare output directory
    sub_str = sub_num2str(sub);
    save_path = fullfile(par.deriv_path, 'erp', sub_str);

    if ~exist(save_path)
        mkdir(save_path);
    end


    %% load and preprocess
    % load data after artefact rejection
    fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub));
    fname =  sprintf('sub-%03d_task-EEG_filt_ep_AB-100_chanint', sub);
    [header, data] = CLW_load(fullfile(fpath, fname));


    % low-pass filter for ERPs
    [header, data] = RLW_butterworth_filter(header, data, ...
        'filter_type', 'lowpass', ...
        'high_cutoff', par.filt_erp_cut_low, ...
        'filter_order', 4);


    % re-segment
    codes = unique({header.events.code});
    [header, data] = segment_safe(header, data, codes,...
        'x_start', 0, 'x_duration', par.trial_dur);


    % rereference to average mastoids
    [header, data] = RLW_rereference(header, data, ...
        'apply_list', {header.chanlocs.labels}, ...
        'reference_list', par.eeg_ref_chans);


    % keep only frontal channels
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
    else
        error('roi name not recognised, ending the script');
    end

    [header, data] = RLW_arrange_channels(header, data, roi_chans);


    % get erp chunks
    [header_chunk, data_chunk] = chunk_cycles_per_cond(...
        header, data, par.cycle_dur, par.experiment);


    % downsample data to 100 Hz (factor 10)
    [header_chunk, data_chunk] = RLW_downsample(header_chunk, data_chunk, ...
        'x_downsample_ratio', 10);
    if (1/header_chunk.xstep) ~= 100
        error('error: check sampling rate and adjust downsampling accordingly');
    end


    % baseline correction (subtract mean of each epoch)
    data_chunk = data_chunk - mean(data_chunk, 6);
    %     size(data_chunk)


    %% chunk rejection based on threshold
    % a chunk is excluded id a range of threshold_range_chunk is exceeded
    % in at least 25% of the roi channels

    chunk_mask = ones(size(data_chunk, 1), 1);

    for i_chunk = 1:size(data_chunk, 1)

        for i_channel = 1:size(data_chunk, 2)
            range(i_channel) = max(data_chunk(i_chunk,i_channel,:,:,:,:))-...
                min(data_chunk(i_chunk, i_channel,:,:,:,:));
        end

        if sum(range>threshold_range_chunk) >= floor(length(roi_chans)*0.25)
            chunk_mask(i_chunk) = 0;
        end
    end

    % remove chunks based on the mask
    [header_chunk, data_chunk] = RLW_arrange_epochs( ...
        header_chunk, data_chunk, find(chunk_mask));
    %     size(data_chunk)


    %% (optional) equalise number chunks by condition

    % count number of chunks per condition
    for i_cond = 1:par.n_cond
        chunk_per_cond(i_sub, i_cond) = sum(strcmp( ...
            {header_chunk.events.code}, num2str(i_cond)));

    end

    % equalise number of trials per condition at the individual level
    num_chunk_equalise = min(chunk_per_cond(i_sub, :));

    epoch_keep = [];
    rng(i_sub); % seed for reproducibility of chosen epochs 
    
    for i_cond = 1:par.n_cond

        % random numbers from 1 to n -numb. chunks for that conditions
        rand_chunk = randperm(chunk_per_cond(i_sub,i_cond));

        % find first event for trigger 'i_cond'
        first_epoch_cond = min(find(strcmp({header_chunk.events.code}, ...
            num2str(i_cond))));

        epoch_keep = [epoch_keep ...
            rand_chunk(1:num_chunk_equalise)+first_epoch_cond-1];

    end

    epoch_keep = sort(epoch_keep);

    [header_chunk, data_chunk] = RLW_arrange_epochs( ...
        header_chunk, data_chunk, epoch_keep);

    for i_cond = 1:par.n_cond
        chunk_per_cond_sanity_check(i_sub, i_cond) = sum(strcmp( ...
            {header_chunk.events.code}, num2str(i_cond)));

    end
    disp(sprintf('%s', sub_num2str(sub)));
    disp(chunk_per_cond_sanity_check);


    % average cycle chunks separately for each condition
    [header_erp_avg, data_erp_avg] = average_trials_per_cond( ...
        header_chunk, data_chunk, [], 'concatenate_epochs', false);


    % save erp data
    fname = sprintf('sub-%03d_response-erp_lpf-%.0fHz_rr-m', ...
        sub, par.filt_erp_cut_low);

    header_erp_avg.name = fname;
    CLW_save(save_path, header_erp_avg, data_erp_avg);


    % sanity check plot
    plot_erp_overlay(squeeze(mean(data_erp_avg, 2)), 1/header_erp_avg.xstep);
    title(sub_num2str(sub));

end