% script to perform grand average

clc
clear all
close all


%% parameters
par = get_par('experiment', 'newborn'); %'newborn' '2ioi_conga'

switch par.experiment

    case 'newborn'

        subjects = [1 5 7 8 9 10 12 17 18 22 23 27 28 31 33 35 36 37 38];
        subj_grand_avg_suffix = sprintf('grand-%d',length(subjects));

    case '2ioi_conga'

        subjects = [2:7 9:15 17:20];
        subj_grand_avg_suffix = 'grand';

end

signal_to_grand_average = 'erp'; %'lf'; %



%%
if strcmp(signal_to_grand_average, 'lf')


    data_lf_all = cell(1, length(subjects));

    for i_sub=1:length(subjects)

        sub = subjects(i_sub);

        sub_str = sub_num2str(sub)


        % load data
        fpath = fullfile(par.deriv_path, 'preproc', sub_str);
        if strcmp(par.experiment, 'newborn')
            fname =  sprintf('sub-%03d_task-EEG_filt_ep_AB-100_chanint_rr-m.lw6', sub);
        else
            fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*_rr-m.lw6', sub);
        end
        d = dir(fullfile(fpath, fname));

        [header, data] = CLW_load(fullfile(d.folder, d.name));

        % save individual files in structure for grand average later
        data_lf_all{i_sub} = data;
    end

    %%
    header_ga = header;
    data_ga = mean(cat(7, data_lf_all{:}), 7);

    save_path = fullfile(par.deriv_path, 'preproc', ...
        sprintf('sub-%s', subj_grand_avg_suffix));

    if ~exist(save_path)
        mkdir(save_path);
    end

    if strcmp(par.experiment, 'newborn')
        header_ga.name = sprintf( ...
            'sub-%s_task-EEG_filt_ep_AB-100_chanint_rr-m', ...
            subj_grand_avg_suffix);
    else
        header_ga.name = sprintf( ...
            'sub-%s_task-EEG_filt_ep_chanint_icfilt_rr-m', ...
            subj_grand_avg_suffix);
    end

    CLW_save(save_path, header_ga, data_ga);


elseif strcmp(signal_to_grand_average, 'erp')


    header_erp_all = cell(1, length(subjects));
    data_erp_all = cell(1, length(subjects));


    for i_sub=1:length(subjects)

        sub = subjects(i_sub);
        sub_str = sub_num2str(sub)

        % load data
        fpath = fullfile(par.deriv_path, 'erp', sub_str);
        fname =  sprintf('sub-%03d_response-erp_lpf-%.0fHz_rr-m.lw6', sub, par.filt_erp_cut_low);

        d = dir(fullfile(fpath, fname));

        [header_erp_ind, data_erp_ind] = CLW_load(fullfile(d.folder, d.name));

        [header_erp_ind.events.sub] = deal(sub);

        header_erp_all{i_sub} = header_erp_ind;
        data_erp_all{i_sub} = data_erp_ind;

    end


    data_erp = cat(1, data_erp_all{:});

    events = cellfun(@(x) x.events, header_erp_all, 'uni', 0);
    header_erp = header_erp_all{1};
    header_erp.events = cat(2, events{:});
    header_erp.datasize = size(data_erp);

    [header_erp_mean, data_erp_mean] = average_trials_per_cond(...
        header_erp, data_erp, [], ...
        'concatenate_epochs', false);

    save_path = fullfile(par.deriv_path, 'erp', ...
        sprintf('sub-%s', subj_grand_avg_suffix));

    if ~exist(save_path)
        mkdir(save_path);
    end

    fname =  sprintf('sub-%s_response-erp_lpf-%.0fHz_rr-m', ...
        subj_grand_avg_suffix, par.filt_erp_cut_low);
    header_erp_mean.name = fname; 
    CLW_save(save_path, header_erp_mean, data_erp_mean); 

    fname =  sprintf('sub-%s_response-erp_lpf-%.0fHz_rr-m_ind', ...
        subj_grand_avg_suffix, par.filt_erp_cut_low);
    header_erp.name = fname;
    CLW_save(save_path, header_erp, data_erp);

    %sanity check (pretty)
    load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'))
    plot_erp_overlay(squeeze(mean(data_erp_mean, 2)), 1/header_erp.xstep, ...
        'ratios', stim.ratios)

end


%% save participants list

participant_list_filename = sprintf('sub-%s_participants_list.txt', ...
    subj_grand_avg_suffix);

fid = fopen(fullfile(save_path, participant_list_filename), 'w');  % Open file for writing

fprintf(fid, 'participants\n');      % Write header

fprintf(fid, '%d\n', subjects);             % Write each element on a new line

fclose(fid);