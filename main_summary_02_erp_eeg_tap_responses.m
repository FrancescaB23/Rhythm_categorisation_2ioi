% script to perform group level statistic and summary of the eeg rsa
% results, erp responses

clear
close all
clc


%% parameters
par = get_par('experiment', '2ioi_conga');
subjects = [num2cell([2:7 9:15 17:20])];

response = 'tap-force-erp';%'erp'; %

roi_name = 'front';
avg_chans = true;

feat_name_stim = 'ratio';

lpf_cutoff = 10;

run_stats = false;


%% prepare file and folder names dependent on response type, output path

switch response
    case 'erp'

        subdirs = fullfile(sprintf('response-%s_stim-%s', response, feat_name_stim), ...
            sprintf('lpf-%.0fHz', lpf_cutoff), ...
            sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)));

        save_path = fullfile(par.summary_path, subdirs);
        mkdir(save_path);

        fname_save = fullfile(save_path, ...
            sprintf('response-%s_lpf-%.0fHz_roi-%s_avgChans-%s_cat-2-3-4-5', ...
            response, ...
            lpf_cutoff, ...
            roi_name, ...
            jsonencode(avg_chans)));

        fname_tbl = sprintf('response-%s_lpf-%.0fHz_roi-%s_avgChans-%s_cat-2-3-4-5.csv', ...
            response, ...
            lpf_cutoff, ...
            roi_name, ...
            jsonencode(avg_chans));

    case 'tap-force-erp'

        subdirs = fullfile(sprintf('response-%s_stim-%s', response, feat_name_stim), ...
            sprintf('lpf-%.0fHz', lpf_cutoff));

        save_path = fullfile(par.summary_path, subdirs);
        mkdir(save_path);

        fname_save = fullfile(save_path, ...
            sprintf('response-%s_lpf-%.0fHz_cat-2-3-4-5', ...
            response, ...
            lpf_cutoff));

        fname_tbl = sprintf('response-%s_lpf-%.0fHz_cat-2-3-4-5.csv', ...
            response, ...
            lpf_cutoff);

end


%% load data

tbl = readtable(fullfile(par.results_path, subdirs, fname_tbl));

res = [];

for i_sub=1:length(subjects)

    res(i_sub).sub = subjects{i_sub};
    sub_str = sub_num2str(res(i_sub).sub);
    load_path = fullfile(par.results_path, subdirs, sub_str);

    fname = fullfile(load_path, ...
        sprintf('%s_feat-erp_method-pearson_rdm.mat', sub_str));

    fname
    
    rdm = load(fname);
    res(i_sub).rdm = rdm.rdm;

    fname = fullfile(load_path, ...
        sprintf('%s_feat-erp_method-pearson_nCateg-2-3-4-5_perm.mat', sub_str));

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

end

% save data
save([fname_save, '_data.mat'], 'res');


%% load stimulus

load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'));

rdm_stim = rsa_stim_lf(par, [], 'feat', feat_name_stim);


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
%             ax.XTick = [1 : 33];
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


% distribution of r values across all tested models
ax = pnl(2, 1, 2, 3).select();
plot_r_over_all_models([res.r_all_models]', [1:33], 'ax', ax);


% group-level permutation test
ax = pnl(2, 1, 1, 2).select();
h = histogram(ax, res_group.r_null_dist, 'DisplayStyle', 'stairs', 'linew', 2);
hold(ax, 'on');
plot(ax, [res_group.r_obs, res_group.r_obs], [0, max(h.Values)], ':r', 'linew', 2)
title(sprintf('r = %.2f, p = %.1g', res_group.r_obs, res_group.p_perm));
xlabel('r');


% grand average RDM
load_path = fullfile(par.results_path, subdirs, 'sub-grand');
fname = 'sub-grand_feat-erp_method-pearson_rdm.mat';

res_grand = load(fullfile(load_path, fname));

ax = pnl(2, 1, 1, 1).select();
plot_rdm(res_grand.rdm, 'ax', ax);
pnl(2, 1, 1, 1).title('grand avg');


% plot individual RDMs
pnl_ind_rdm = pnl(1);

plot_all_rdms(pnl_ind_rdm, ...
    [res.sub], ...
    {res.rdm}, ...
    {res.best_model}, ...
    [res.p_perm] ...
    );

pnl_ind_rdm.de.margin = [0, 0, 0, 8];
pnl_ind_rdm.margin = [20, 2, 2, 15];

pnl.marginleft = 15;
pnl.margintop = 15;
pnl.fontsize = 12;


print(f, '-dpng', '-painters', '-r600', [fname_save, '_summary.png']);



%% save summary report text file

fid = fopen([fname_save, '_readme.txt'], 'w');

fprintf(fid, 'mean r = %.4f (SD = %.4f)\n\n', mean([res.r]), std([res.r]));

fprintf(fid, 'permuted pvalue (group level) = %.3g \n\n', res_group.p_perm);

fprintf(fid, 'categorical model significant in %d subjects (out of %d): \n', ...
    sum([res.p_perm] < 0.05), length(res));

fprintf(fid, '%s \n\n', ...
    vec_to_str([subjects{[res.p_perm] < 0.05}], 'format','%d', 'sep','\n'));

fclose(fid);


%%
close all