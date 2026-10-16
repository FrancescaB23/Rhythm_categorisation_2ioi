% compare RSA outputs obtain from tap magnitude spectra and ITI ratio
% responses

clear
close all
clc


%% parameters

par = get_par('experiment', '2ioi_conga');
subjects = [num2cell([2:7 9:15 17:20])];

feat_name_stim = 'ratio';

feat_snr = false;
harmonics = [1:6];

n_perm_rsa_group = 10000;

run_stats = false;
fig_save = false;


%%
if feat_snr
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2));
else
    snr_str = 'none';
end

harm_str = vec_to_str(harmonics, 'format','%d', 'sep','-');

save_path = fullfile(par.deriv_path, 'comparisons', ...
    'tap_mX_vs_ITIratio',...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str));
if ~exist(save_path)
    mkdir(save_path);
end


%% load data for mX
response = 'tap-force';
feat_name_resp = 'mX';

subdirs = fullfile(sprintf('response-%s_stim-%s', response, feat_name_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str));

fname = fullfile(par.summary_path, subdirs, ...
    sprintf('response-%s_feat-%s_harm-%s_snr-%s_cat-2-3-4-5_data.mat', ...
    response, ...
    feat_name_resp, ...
    harm_str, ...
    snr_str));

res_mX = load(fname);


%% load data for ITI ratio
response = 'tap-onsets';

subdirs = fullfile(sprintf('response-%s_stim-%s', ...
    response, feat_name_stim));

fname = fullfile(par.summary_path, subdirs, ...
    sprintf('response-%s_stim-%s_data.mat', ...
    response, ...
    feat_name_stim));

res_ratio = load(fname);


%% sanity check
if isequal([res_mX.res.sub], [res_ratio.res.sub], [subjects{:}])
    disp('subject order across the two measures match, continuing the script');
else
    error('subject order across the two measures do not match, check!');
end


%% save fname
fname_save_base = fullfile(save_path, ...
    sprintf('response-tap-mX-versus-ratio_featStim-%s_harm-%s_snr-%s_feat-mXvsRatio', ...
    feat_name_stim, ...
    harm_str, ...
    snr_str));


%% correlate RDMs based on mX and tap ratio

% load stimulus
load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'));
rdm_stim = rsa_stim_lf(par, [], 'feat', feat_name_stim);

% group-level permutation test
if exist([fname_save_base, '_permGroup.mat']) && ~run_stats
    clear res_group
    load([fname_save_base, '_permGroup.mat']);
else
    res_group = [];
    [res_group.r_obs, res_group.r_null_dist, res_group.p_perm] = ...
        perm_rdm_corr_group_empirical_model(...
        {res_mX.res.rdm}, ...
        {res_ratio.res.rdm},  ...
        n_perm_rsa_group, ...
        'parfor', true, ...
        'covariate', rdm_stim);
    save([fname_save_base, '_permGroup.mat'], 'res_group');
end

% individual level permutation test with fixed model
if exist([fname_save_base, '_permInd.mat']) && ~run_stats
    clear res_ind
    load([fname_save_base, '_permInd.mat']);
else
    res_ind = [];
    for i_sub=1:length(subjects)
        i_sub
        [res_ind(i_sub).r_obs, res_ind(i_sub).p_perm, res_ind(i_sub).r_null_dist] = ...
            perm_rdm_corr(...
            res_mX.res(i_sub).rdm, ...
            res_ratio.res(i_sub).rdm, ...
            n_perm_rsa_group, ...
            'parfor', true, ...
            'covariate', rdm_stim);
    end
    save([fname_save_base, '_permInd.mat'], 'res_ind');
end

% RHO of individual subjects
f = figure('pos', [809 696 80 173], 'color', 'white');
pnl = panel(f);
ax = pnl.select();

plot_r_bar([res_ind.r_obs], [res_ind.p_perm], 'y_maxlim', 1, 'ax', ax);

ax.YLim = [-1, 1];
ax.YTick = [-1, 1];
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
txt = text(0, 1, sprintf('p=%.1g', res_group.p_perm));
txt.HorizontalAlignment = 'center';
pnl.margin = [8,4,0,1];

disp(sprintf(['\n\npartial correlation tap mX-ratio, rho = %.4f, p = %.4f\n' ...
    'significant in %d participants out of %d (p<0.05)'], ...
    res_group.r_obs, res_group.p_perm, ...
    sum([res_ind.p_perm]<.05), length(subjects)));


if fig_save
    saveas(f, [fname_save_base, '_corrRdmsIndRho.png']);
end


%% compare rho with best model

r_mX = [res_mX.res.r]';
r_ratio = [res_ratio.res.r]';

z_mX = r_to_z(r_mX);
z_ratio = r_to_z(r_ratio);

c = [135, 175, 201]/255;

f = figure('pos', [1517 809 144 193], 'color', 'white');
pnl = panel(f);
ax = pnl.select();

s = scatter(repmat(0, 1, length(z_mX)), z_mX, 60, c, 'filled');
hold on ;
s = scatter(repmat(1, 1, length(z_ratio)), z_ratio, 60, c, 'filled');
plot(repmat([0;1],1,length(z_mX)), [z_mX, z_ratio]', '-', 'linew', 2, 'color', c);
ax.XLim = [-0.5, 1.5];

ax.XTick = [0, 1];
ax.XTickLabel = {'mX', 'ratio'};
pnl.margin = [10,10,10,10];

[h, p, ci, stats] = ttest(z_mX, z_ratio);

bf10 = bf.ttest(z_mX, z_ratio);   % paired samples BF

txt = text(0.5, 1.1*max([z_mX; z_ratio]), ...
    sprintf('r_mX=%.4f r_ratio=%.4f\nBF10=%.4f, \np=%.4g', ...
    mean(r_mX), mean(r_ratio), bf10, p));
txt.HorizontalAlignment = 'center';

disp(sprintf(['\n\nttest mX-ratio (Fisher transformed), ' ...
    'p = %.4f\nBF10 = %.4f'], p, bf10));

if fig_save
    saveas(f, sprintf('%s_ttestRhoFisher.png', fname_save_base));
end


%% distribution of correlations across categorical model boundaries

r_mX_distr = [res_mX.res.r_all_models]';
r_ratio_distr = [res_ratio.res.r_all_models]';

z_mX_distr = r_to_z(r_mX_distr);
z_ratio_distr = r_to_z(r_ratio_distr);

r = diag(corr(z_mX_distr', z_ratio_distr', 'type', 'pearson'));
z = r_to_z(r);

c_group_mX = [31, 119, 180] / 255;      % Blue (for mean)
c_ind_mX = [174, 199, 232] / 255;    % Light blue (for std)

c_group_ratio = [255, 127, 14] / 255;      % Orange (for mean)
c_ind_ratio = [255, 187, 120] / 255;    % Light orange (for std)

f = figure('pos', [476 675 358 191], 'color', 'white');
pnl = panel(f);
pnl.pack('h', [86, 14]);

ax = pnl(1).select();

tested_bounds = res_mX.res(1).all_models_info.mod_id';

mu_mX = mean(r_mX_distr);
sem_mX = std(r_mX_distr, [], 1) ./ sqrt(size(r_mX_distr, 1));

mu_ratio = mean(r_ratio_distr);
sem_ratio = std(r_ratio_distr, [], 1) ./ sqrt(size(r_ratio_distr, 1));

x = [tested_bounds, fliplr(tested_bounds)];
y = [mu_ratio + sem_ratio, fliplr(mu_ratio - sem_ratio)];
fill(x, y, c_ind_ratio, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
plot(tested_bounds, mu_ratio, 'LineWidth', 3.5, 'Color', c_group_ratio);
hold on;

x = [tested_bounds, fliplr(tested_bounds)];
y = [mu_mX + sem_mX, fliplr(mu_mX - sem_mX)];
fill(x, y, c_ind_mX, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
plot(tested_bounds, mu_mX, 'LineWidth', 3.5, 'Color', c_group_mX);
ax.XAxis.TickLength = [0,0];
ax.YAxis.TickLength = [0,0];
ax.XTick = [];
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ylabel('r');

xlim([.5 33.5]);
ylim([-0.7 0.75]);


ax = pnl(2).select();

c = [135, 175, 201] / 255;
plot(0 + 0.2*rand(1,length(z))-0.1, z, 'o', 'color', c, 'markerfacecolor', c, 'markersize', 5);

z_mu = mean(z);
z_sem = std(z) / sqrt(length(z));
z_ci = norminv(1 - 0.025) * z_sem;

hold on
c = [0 0 0];
plot(0, z_mu, 'o', 'color', c, 'markerfacecolor', c, 'markersize', 8);
plot([0,0], [z_mu-z_ci, z_mu+z_ci], 'color', c, 'linew', 3);

ax.XLim = [-0.2, 0.2];
ax.XTick = [];
ax.XAxisLocation = 'origin';
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.XAxis.TickLength = [0,0];
ax.YAxis.TickLength = [0,0];
ax.YLim = [-1, 3];
% ax.YTick = sort([min(ax.YTick), 0, 4]);
ylabel('Fisher-transformed r');

[h, p, ci, stats] = ttest(z);

txt = text(0, max(z)*1.2, sprintf('r = %.2f\np = %.2g', mean(r), p));
txt.HorizontalAlignment = 'center';

if fig_save
    saveas(f, [fname_save_base, '_distrCorr.png']);
end


disp(sprintf('\n\ncorrelation distribution rho mX-ratio, r = %.4f p = %.6f', ...
    mean(r), p));


%% figures for manuscript
% correlation coefficients and overlay best models for tap ratio

fig_path = fullfile(par.deriv_path, 'figures', 'tap_ratio');
if ~exist(fig_path)
    mkdir(fig_path);
end

c_group = [136, 34, 85] / 255;
c_ind = [204, 102, 119] / 255;

f = figure('pos', [827 508 153 189], 'color', 'white');

pnl = panel(f);
ax = pnl.select();

plot_r_bar([res_ratio.res.r], [res_ratio.res.p_perm], ...
    'ax', ax, 'c_group', c_group, 'c_ind', c_ind);

ax.YLim = [0, 1];
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.YTick = [0, 1];
ax.FontSize = 14; 
ax.FontName = 'Arial';
ylabel(ax, 'ρ', 'FontName', 'Arial', 'FontSize', 14);

if fig_save
    fname_save = sprintf('rho-adults_TAP-ratio_stim-%s', ...
        feat_name_stim);
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_path, fname_save), '.png']);
end


% plot overlay best categorical model, all participants
f = plot_rdms_overlay({res_ratio.res.best_model}, 'alpha', 0.1);
title('TAP ratio all sub');

if fig_save
    fname_save = sprintf('rdm-adults-best-model_TAP-ratio_stim-%s', ...
        feat_name_stim);
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_path, fname_save), '.png']);
end

