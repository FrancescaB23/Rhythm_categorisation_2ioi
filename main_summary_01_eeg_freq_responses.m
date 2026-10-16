% script to perform group level statistic and summary of the eeg rsa
% results, frequency domain

clear
close all
clc


%% parameters
par = get_par('experiment', '2ioi_conga'); %'newborn' '2ioi_conga'

feat_names =  {'mX'};%{'X', 'mX'}; %

roi_name = 'front';
avg_chans = true;

feat_snr = false;

harmonics = [1:6];

% !!!!!!!!!!!! NOTE
% Check feat_name_stim in the script

run_stats = false;


%%

if strcmp(par.experiment, 'newborn')

    subjects = num2cell([1 5 7 8 9 10 12 17 18 22 23 27 28 31 33 35 36 37 38]);
    ga_code = sprintf('%d', length(subjects));

elseif strcmp(par.experiment, '2ioi_conga')

    subjects = [num2cell([2:7 9:15 17:20])];
end


%%
for i_feat=1:length(feat_names)


    feat_name_resp = feat_names{i_feat};

%     feat_name_stim = feat_names{i_feat};
    feat_name_stim = 'ratio';
%     feat_name_stim = 'cochlear';


    %% prepare output paths
    harm_str = vec_to_str(harmonics, 'format','%d', 'sep','-');

    if feat_snr
        snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2));
    else
        snr_str = 'none';
    end

    subdirs = fullfile(sprintf('response-lf_stim-%s', feat_name_stim), ...
        sprintf('harm-%s', harm_str), ...
        sprintf('snr-%s', snr_str), ...
        sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)) ...
        );

    save_path = fullfile(par.summary_path, subdirs);
    mkdir(save_path);

    fname_save = fullfile(save_path, ...
        sprintf('response-lf_feat-%s_harm-%s_snr-%s_roi-%s_avgChans-%s_cat-2-3-4-5', ...
        feat_name_resp, ...
        harm_str, ...
        snr_str, ...
        roi_name, ...
        jsonencode(avg_chans) ...
        ));


    %% load data

    fname_tbl = sprintf('response-lf_feat-%s_harm-%s_snr-%s_roi-%s_avgChans-%s_cat-2-3-4-5.csv', ...
        feat_name_resp, ...
        harm_str, ...
        snr_str, ...
        roi_name, ...
        jsonencode(avg_chans) ...
        );

    tbl = readtable(fullfile(par.results_path, subdirs, fname_tbl));

    res = [];

    for i_sub=1:length(subjects)

        res(i_sub).sub = subjects{i_sub};
        sub_str = sub_num2str(res(i_sub).sub);

        load_path_sub = fullfile(par.results_path, subdirs, sub_str);

        fname = fullfile(load_path_sub, ...
            sprintf('%s_feat-%s_method-pearson_rdm.mat', sub_str, feat_name_resp));

        rdm = load(fname);
        res(i_sub).rdm = rdm.rdm;

        fname = fullfile(load_path_sub, ...
            sprintf('%s_feat-%s_method-pearson_nCateg-2-3-4-5_perm.mat', sub_str, feat_name_resp));

        perm_res = load(fname);

        res(i_sub).r = perm_res.r_obs;
        res(i_sub).p_perm = perm_res.p_perm;
        res(i_sub).best_model = perm_res.best_model;

        [res(i_sub).all_models_info, idx] = ...
            sortrows(perm_res.all_models_info, ...
            {'n_categ', 'bound1', 'bound2', 'bound3', 'bound4'});
        res(i_sub).all_models_info.mod_id = [1 : length(res(i_sub).all_models_info.n_categ)]';
        res(i_sub).r_all_models = res(i_sub).all_models_info.r;
        [tmp_max_r, res(i_sub).best_model_id] = max(res(i_sub).all_models_info.r);
        if ~isequal(tmp_max_r, res(i_sub).r)
            disp('check code - max indr does not correspond')
        end
        res(i_sub).r_null_dist = perm_res.r_null_dist;
        mask = tbl.sub == res(i_sub).sub;
        res(i_sub).z_snr = tbl{mask, 'z_snr'};

    end

    % save data
    save([fname_save, '_data.mat'], 'res');


    %% load stimulus

    load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'));

    if strcmp(feat_name_stim, 'cochlear') % cochlear model as regressor

        % load cochlear fft responses
        fpath_cochlear = fullfile(['/Users/barbero/Documents/Projects/' ...
            'Rhototype_2ioi_newborns/adults/data/derivatives'], ...
            'fft', 'response-urear-an');
        fname_cochlear = sprintf('response-urear-an_snr-none_%s', feat_name_resp);
        [header_cochlear, data_cochlear] = CLW_load(fullfile(fpath_cochlear, fname_cochlear));

        % get frequency vector - specific for cochlear responses
        freq_cochlear = [0: header_cochlear.datasize(end)-1] * header_cochlear.xstep;
        frex = 1 / par.cycle_dur * harmonics; 

        % get cochlear rdm
        feat_for_rdm_cochlear = get_feat_from_fft(header_cochlear, data_cochlear, freq_cochlear, frex);

        rdm_stim = get_rdm_from_feat(feat_for_rdm_cochlear, 'method', 'pearson');

        figure()
        imagesc(rdm_stim)
        title(sprintf('cochlear %s', feat_name_resp))

    else

        rdm_stim = rsa_stim_lf(par, 1/par.cycle_dur * harmonics, 'feat', feat_name_stim);

        figure()
        imagesc(rdm_stim);
        title(sprintf('%s', feat_name_stim));
    end


    %% group-level permutation test
    if exist([fname_save, '_permGroup.mat']) && ~run_stats
        clear res_group
        load([fname_save, '_permGroup.mat']);
    else
        res_group = [];

        [res_group.r_obs, res_group.r_null_dist, res_group.p_perm] = ...
            perm_rdm_corr_group(...
            {res.rdm}, [2, 3, 4, 5], 10000, ...
            'covariate', rdm_stim, ...
            'n_mini_categ_skip', 1);

        save([fname_save, '_permGroup.mat'], 'res_group');
    end


    %% plotting

    f = figure('pos', [272 1 984 1046], 'color', 'white');

    pnl = panel(f);

    pnl.pack('h', [20, 50, 30]);

    pnl(2).pack('v', [25, 25, 50]);

    pnl(2, 1).pack('v', [50, 50]);
    pnl(2, 1, 1).pack('h', 3);
    pnl(2, 1, 2).pack('h', 3);

    pnl(2, 2).pack('v', 3);


    % bar plot with individual subjects overlaid as points
    ax = pnl(2, 1, 2, 1).select();
    plot_r_bar([res.r], [res.p_perm], 'y_maxlim', 1, 'ax', ax);
    ax.YLabel.String = 'r';

    % bar plot with zSNR of individual subjects overlaid as points
    ax = pnl(2, 1, 2, 2).select();
    plot_r_bar([res.z_snr], [res.p_perm], 'y_maxlim', Inf, 'ax', ax);
    ax.YTick = [0, round(mean([res.z_snr]), 1), floor(max([res.z_snr])*10)/10];
    ax.YLabel.String = 'z snr';

    % best fitting models across all subjects
    pnl(2,2,1).pack('h', [30, 70]);
    ax = pnl(2, 2, 1, 1).select();
    plot_rdms_overlay({res.best_model}, 'alpha', 0.1, 'ax', ax);
    ax = pnl(2, 2, 1, 2).select();
    col = [0.6, 0.6, 0.6];
    histogram(ax, [res.best_model_id], ...
        [0.5 : 1 : 33+0.5], ...
        'EdgeColor', col, 'FaceColor', col, 'LineWidth', 2);
    box(ax, 'off');
    ax.YTick = [0 : 1 : ceil(ax.YLim(2))];
    ax.FontSize = 12;
    pnl(2, 2, 1).title('all');

    % best fitting models across significant subjects
    idx_signif = [res.p_perm] < 0.05;
    pnl(2,2,2).pack('h', [30, 70]);
    ax = pnl(2, 2, 2, 1).select();
    plot_rdms_overlay({res(idx_signif).best_model}, 'alpha', 0.1, 'ax', ax);
    ax = pnl(2, 2, 2, 2).select();
    col = [0.6, 0.6, 0.6];
    histogram(ax, [res(idx_signif).best_model_id], ...
        [0.5 : 1 : 33+0.5], ...
        'EdgeColor', col, 'FaceColor', col, 'LineWidth', 2);
    box(ax, 'off');
    ax.YTick = [0 : 1 : ceil(ax.YLim(2))];
    ax.FontSize = 12;
    pnl(2, 2, 2).title('signif categ');

    % best fitting models across subjects with significant z-snr
    idx_signif = [res.z_snr] > norminv(1 - 1e-3);
    pnl(2,2,3).pack('h', [30, 70]);
    ax = pnl(2, 2, 3, 1).select();
    plot_rdms_overlay({res(idx_signif).best_model}, 'alpha', 0.1, 'ax', ax);
    ax = pnl(2, 2, 3, 2).select();
    col = [0.6, 0.6, 0.6];
    histogram(ax, [res(idx_signif).best_model_id], ...
        [0.5 : 1 : 33+0.5], ...
        'EdgeColor', col, 'FaceColor', col, 'LineWidth', 2);
    box(ax, 'off');
    ax.YTick = [0 : 1 : ceil(ax.YLim(2))];
    ax.FontSize = 12;
    pnl(2, 2, 3).title('signif zsnr');
    pnl(2, 2, 3, 2).xlabel('boundary location');


    % distribution of r values across all tested models
    ax = pnl(2, 1, 2, 3).select();
    plot_r_over_all_models([res.r_all_models]', [1:33], 'ax', ax);


    % relationhsip between z-snr and r
    ax = pnl(2, 1, 1, 3).select();
    plot(ax, [res.z_snr], [res.r], 'ko')
    axis(ax, 'square');
    box(ax, 'off');
    ax.XLabel.String = 'z snr';
    ax.YLabel.String = 'r';


    % group-level permutation test
    ax = pnl(2, 1, 1, 2).select();
    h = histogram(ax, res_group.r_null_dist, 'DisplayStyle', 'stairs', 'linew', 2);
    hold(ax, 'on');
    plot(ax, [res_group.r_obs, res_group.r_obs], [0, max(h.Values)], ':r', 'linew', 2)
    title(sprintf('r = %.2f, p = %.1g', res_group.r_obs, res_group.p_perm));
    xlabel('r');


    % grand average RDM
    if strcmp(par.experiment, 'newborn')
        load_path = fullfile(par.results_path, subdirs, ...
            sprintf('sub-grand-%s', ga_code));
        fname = sprintf('sub-grand-%s_feat-%s_method-pearson_rdm.mat', ...
            ga_code, feat_name_resp);
    else
        load_path = fullfile(par.results_path, subdirs, 'sub-grand');
        fname = sprintf('sub-grand_feat-%s_method-pearson_rdm.mat', feat_name_resp);
    end

    res_grand = load(fullfile(load_path, fname));

    ax = pnl(2, 1, 1, 1).select();
    plot_rdm(res_grand.rdm, 'ax', ax);
    if strcmp(par.experiment, 'newborn')
        pnl(2, 1, 1, 1).title(sprintf('grand avg %s', ga_code));
    else
        pnl(2, 1, 1, 1).title('grand avg');
    end


    % grand average mX
    if strcmp(par.experiment, 'newborn')
        load_path = fullfile(par.results_path, subdirs, ...
            sprintf('sub-grand-%s', ga_code));
        fname = sprintf('sub-grand-%s_mX.mat', ga_code);
    else
        load_path = fullfile(par.results_path, subdirs, 'sub-grand');
        fname = 'sub-grand_mX.mat';
    end
    mX_data = load(fullfile(load_path, fname));

    pnl_mX = pnl(3);

    plot_mX_lw(mX_data.header_mX_to_plot, mX_data.data_mX_to_plot, ...
        [1 : mX_data.header_mX_to_plot.datasize(2)], ...
        mX_data.frex, ...
        'z_snr', mX_data.z_snr, ...
        'fmax', max([mX_data.frex, 15]), ...
        'pnl', pnl_mX);


    % plot individual RDMs
    pnl_ind_rdm = pnl(1);


    % if many subjects, all individual data won't fit in the panel,
    % just plot a subset of those
    if length(res) > 20

        res_for_panel = res(1:20);
        plot_all_rdms(pnl_ind_rdm, ...
            [res_for_panel.sub], ...
            {res_for_panel.rdm}, ...
            {res_for_panel.best_model}, ...
            [res_for_panel.p_perm], ...
            [res_for_panel.z_snr] ...
            )

    else
        plot_all_rdms(pnl_ind_rdm, ...
            [res.sub], ...
            {res.rdm}, ...
            {res.best_model}, ...
            [res.p_perm], ...
            [res.z_snr] ...
            )
    end


    pnl_ind_rdm.de.margin = [0, 0, 0, 8];
    pnl_ind_rdm.margin = [20, 2, 2, 15];

    pnl.marginleft = 15;
    pnl.margintop = 15;
    pnl.fontsize = 12;



    print(f, '-dpng', '-painters', '-r600', [fname_save, '_summary.png']);



    %% save summary report text file

    fid = fopen([fname_save, '_readme.txt'], 'w');

    fprintf(fid, 'mean r = %.4f (SD = %.4f)\n\n', mean([res.r]), std([res.r]));

    fprintf(fid, 'mean zSNR = %.4f (SD = %.4f)\n\n', mean([res.z_snr]), std([res.z_snr]));

    fprintf(fid, 'permuted pvalue (group level) = %.3g \n\n', res_group.p_perm);

    fprintf(fid, 'categorical model significant in %d subjects (out of %d): \n', ...
        sum([res.p_perm] < 0.05), length(res));

    fprintf(fid, '%s \n\n', ...
        vec_to_str([subjects{[res.p_perm] < 0.05}], 'format','%d', 'sep','\n'));

    fclose(fid);


    %%
    close all

end