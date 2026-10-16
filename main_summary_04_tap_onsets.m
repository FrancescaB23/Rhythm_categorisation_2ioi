% script to perform group level statistic and summary of the tap rsa
% results, ITI ratios

clear
close all
clc


%% parameters
par = get_par('experiment', '2ioi_conga');
subjects = num2cell([2:7 9:15 17:20]);

response = 'tap-onsets';
feat_name_stim = 'ratio';

run_stats = false;


%% prepare output path
subdirs = fullfile(sprintf('response-%s_stim-%s', response, feat_name_stim));

save_path = fullfile(par.summary_path, subdirs);
mkdir(save_path);

fname_save = fullfile(save_path, sprintf('response-%s_stim-%s', ...
    response, feat_name_stim));


%% load data

fname_tbl = sprintf('response-%s_cat-2-3-4-5.csv', response);

tbl = readtable(fullfile(par.results_path, subdirs, fname_tbl));

res = [];

for i_sub=1:length(subjects)

    res(i_sub).sub = subjects{i_sub};
    sub_str = sub_num2str(res(i_sub).sub);

    load_path_sub = fullfile(par.results_path, subdirs, sub_str);

    fname = fullfile(load_path_sub, ...
        sprintf('%s_feat-ratios_method-pearson_rdm.mat', sub_str));

    rdm = load(fname);
    res(i_sub).rdm = rdm.rdm;

    fname = fullfile(load_path_sub, ...
        sprintf('%s_feat-ratios_method-pearson_nCateg-2-3-4-5_perm.mat', sub_str));

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

    % let's also load the mean tapped ratio for each subject
    data = load(fullfile(par.deriv_path, 'tap_onsets', sub_str, ...
        sprintf('%s_response-tap-ratios.mat', sub_str)));
    res(i_sub).tap_ratios = data.tap_ratios;

end

% save data
save([fname_save, '_data.mat'], 'res');


%% load stimulus info

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
fname = 'sub-grand_feat-ratios_method-pearson_rdm.mat';

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
    )

pnl_ind_rdm.de.margin = [0, 0, 0, 8];
pnl_ind_rdm.margin = [20, 2, 2, 15];

pnl.marginleft = 15;
pnl.margintop = 15;
pnl.fontsize = 12;

tit = pnl_ind_rdm.title('tap-ratios');
tit.Position(2) = 1.03;
tit.FontWeight = 'bold';

print(f, '-dpng', '-painters', '-r600', [fname_save, '_summary.png']);


%% fit sigmoid to mean tapped ratios across conditions

x = repmat(stim.ratios(:, 1), length(res), 1);
y = cat(1, res.tap_ratios);
y = y(:, 1);

[param, stats, y_pred] = fit_sigmoid(x, y, [], [], 0);

% "min", "max", "x50", "slope"
min_ = param(1);
max_ = param(2);
x50_ = param(3);
slope_ = param(4);

f = figure('color', 'white', 'Position',  [680 797 114 181]);
plot(x, y, 'k.')

hold on
x_vector = [min(x) : (max(x)-min(x))/100 : max(x)];
h = plot(x_vector, y_pred, 'r-', 'linew', 2);
h = plot([x50_, x50_], get(gca, 'YLim'), 'k', 'linew', 4);
h.Color(4) = 0.2;

box off

xlim([min(x)*0.99 max(x)*1.01])

ratio_labs = get_ratio_labs(stim.ratios);

xticks(stim.ratios([1, end], 1));
xticklabels(ratio_labs([1, end]))

yticks(stim.ratios([1, end], 1));
yticklabels(ratio_labs([1, end]))


print(f, '-dpng', '-painters', '-r300', [fname_save, '_sigmFit.png']);



%% save summary report text file

fid = fopen([fname_save, '_readme.txt'], 'w');

fprintf(fid, 'mean r = %.4f (SD = %.4f)\n\n', mean([res.r]), std([res.r]));

fprintf(fid, 'permuted pvalue (group level) = %.3g \n\n', res_group.p_perm);

fprintf(fid, 'categorial model significnat in %d subjects (out of %d): \n', ...
    sum([res.p_perm] < 0.05), length(res));

fprintf(fid, '%s \n\n', ...
    vec_to_str([subjects{[res.p_perm] < 0.05}], 'format','%d', 'sep','\n'));

fclose(fid);


%%
close all;
