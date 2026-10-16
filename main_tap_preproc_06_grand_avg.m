% script to perform grand average tapping responses

clc
clear
close all


%% parameters
par = get_par('experiment', '2ioi_conga');
subjects = [2:7 9:15 17:20];

signal_to_grand_average = 'tap-force-erp'; % 'tap-force';%'ITI_ratio';%


%%
if strcmp(signal_to_grand_average, 'ITI_ratio')

    data_tap_ratios_all = cell(1, length(subjects));


    for i_sub = 1:length(subjects)

        % load data
        sub = subjects(i_sub);

        fpath = fullfile(par.deriv_path, 'tap_onsets', sub_num2str(sub));
        fname =  sprintf('sub-%03d_response-tap-ratios.mat', sub);
        data_tap_ratios = load(fullfile(fpath, fname));

        % store in structure
        data_tap_ratios_all{i_sub} = data_tap_ratios.tap_ratios;
    end


    tap_ratios = mean(cat(3, data_tap_ratios_all{:}), 3);

    save_path = fullfile(par.deriv_path, 'tap_onsets', 'sub-grand');
    mkdir(save_path);

    fname =  sprintf('sub-grand_response-tap-ratios.mat');

    save(fullfile(save_path, fname), 'tap_ratios');


elseif strcmp(signal_to_grand_average, 'tap-force')

    data_tap_all = cell(1, length(subjects));


    for i_sub = 1:length(subjects)

        % load data
        sub = subjects(i_sub);

        fpath = fullfile(par.deriv_path, 'tapping', sub_num2str(sub));
        fname = sprintf('sub-%03d_task-TAP_response-force_avg.lw6', sub);
        d = dir(fullfile(fpath, fname));
        [header_tap, data_tap] = CLW_load(fullfile(d.folder, d.name));

        % store in structure
        data_tap_all{i_sub} = data_tap;
    end

    header_all = header_tap;
    data_all = mean(cat(7, data_tap_all{:}), 7);
    header_all.name = 'sub-grand_task-TAP_response-force_avg';

    save_path = fullfile(par.deriv_path, 'tapping', 'sub-grand');
    mkdir(save_path);

    CLW_save(save_path, header_all, data_all);


elseif strcmp(signal_to_grand_average, 'tap-force-erp')

    header_tap_erp_all = cell(1, length(subjects));
    data_tap_erp_all = cell(1, length(subjects));


    for i_sub=1:length(subjects)

        sub = subjects(i_sub);
        sub_str = sub_num2str(sub)

        % load data
        fpath = fullfile(par.deriv_path, 'erp', sub_str);
        fname =  sprintf('sub-%03d_response-tap-force-erp_lpf-%.0fHz', ...
            sub, par.filt_tap_erp_cut_low);

        [header_tap_erp_ind, data_tap_erp_ind] = ...
            CLW_load(fullfile(fpath, fname));

        [header_tap_erp_ind.events.sub] = deal(sub);

        header_tap_erp_all{i_sub} = header_tap_erp_ind;
        data_tap_erp_all{i_sub} = data_tap_erp_ind;

    end


    data_erp = cat(1, data_tap_erp_all{:});

    events = cellfun(@(x) x.events, header_tap_erp_all, 'uni', 0);
    header_erp = header_tap_erp_all{1};
    header_erp.events = cat(2, events{:});
    header_erp.datasize = size(data_erp);

    [header_tap_erp_mean, data_tap_erp_mean] = average_trials_per_cond(...
        header_erp, data_erp, par, ...
        'concatenate_epochs', false);

    save_path = fullfile(par.deriv_path, 'erp', 'sub-grand');

    if ~exist(save_path)
        mkdir(save_path);
    end

    fname =  sprintf('sub-grand_response-tap-force-erp_lpf-%.0fHz', ...
        par.filt_tap_erp_cut_low);
    header_tap_erp_mean.name = fname;
    CLW_save(save_path, header_tap_erp_mean, data_tap_erp_mean);

    fname =  sprintf('sub-grand_response-tap-force-erp_lpf-%.0fHz_ind', ...
        par.filt_tap_erp_cut_low);
    header_erp.name = fname;
    CLW_save(save_path, header_erp, data_erp);

    %sanity check (pretty)
    load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'))
    plot_erp_overlay(squeeze(mean(data_tap_erp_mean, 2)), 1/header_erp.xstep, ...
        'ratios', stim.ratios)

end
