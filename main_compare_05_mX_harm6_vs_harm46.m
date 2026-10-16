% compare RSA outputs obtained from different harmonics selections

clear
close all
clc


%% parameters
par = get_par('experiment', '2ioi_conga');
subjects = [num2cell([2:7 9:15 17:20])];

measure = 'mX';
harm_selection_A = [1:6];
harm_selection_B = [1:46];

feat_name_stim = 'ratio';
feat_snr = false;
roi_name = 'front';
avg_chans = true;

n_perm_rsa_group = 10000;

run_stats = false;
fig_save = false;


%%
if feat_snr
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2));
else
    snr_str = 'none';
end

harm_str_A = vec_to_str(harm_selection_A, 'format','%d', 'sep','-');
harm_str_B = vec_to_str(harm_selection_B, 'format','%d', 'sep','-');

save_path = fullfile(par.deriv_path, 'comparisons', ...
    sprintf('%s-h%d_vs_h%d', measure, ...
    harm_selection_A(end), harm_selection_B(end)), ...
    sprintf('snr-%s', snr_str));
if ~exist(save_path)
    mkdir(save_path);
end


%% load data for measure A
subdirs = fullfile(sprintf('response-lf_stim-%s', feat_name_stim), ...
    sprintf('harm-%s', harm_str_A), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)));

fname = fullfile(par.summary_path, subdirs, ...
    sprintf('response-lf_feat-%s_harm-%s_snr-%s_roi-%s_avgChans-%s_cat-2-3-4-5_data.mat', ...
    measure, ...
    harm_str_A, ...
    snr_str, ...
    roi_name, ...
    jsonencode(avg_chans)));

res_A = load(fname);


%% load data for measure B
subdirs = fullfile(sprintf('response-lf_stim-%s', feat_name_stim), ...
    sprintf('harm-%s', harm_str_B), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)));

fname = fullfile(par.summary_path, subdirs, ...
    sprintf('response-lf_feat-%s_harm-%s_snr-%s_roi-%s_avgChans-%s_cat-2-3-4-5_data.mat', ...
    measure, ...
    harm_str_B, ...
    snr_str, ...
    roi_name, ...
    jsonencode(avg_chans)));

res_B = load(fname);


%% sanity check
if isequal([res_A.res.sub], [res_B.res.sub], [subjects{:}])
    disp('subject order across the two measures match, continuing the script');
else
    error('subject order across the two measures do not match, check!');
end


%% save fname

fname_save_base = fullfile(save_path, ...
    sprintf('response-%s_featStim-%s_harm-%d-versus-harm-%d_snr-%s_roi-%s_avgChans-%s', ...
    measure, ...
    feat_name_stim, ...
    harm_selection_A(end), ...
    harm_selection_B(end),...
    snr_str, ...
    roi_name, ...
    jsonencode(avg_chans)));


%% correlate RDMs based on measures A and B

% load stimulus
load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'));
rdm_stim = rsa_stim_lf(par,[], 'feat', feat_name_stim);

% group-level permutation test
if exist([fname_save_base, '_permGroup.mat']) && ~run_stats
    clear res_group
    load([fname_save_base, '_permGroup.mat']);
else
    res_group = [];
    [res_group.r_obs, res_group.r_null_dist, res_group.p_perm] = ...
        perm_rdm_corr_group_empirical_model(...
        {res_A.res.rdm}, ...
        {res_B.res.rdm},  ...
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
            res_A.res(i_sub).rdm, ...
            res_B.res(i_sub).rdm, ...
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


disp(sprintf(['\n\npartial correlation %s h-%d - h-%d, rho = %.4f, ' ...
    'p = %.4f\nsignificant in %d participants out of %d (p<0.05)'], ...
    measure, harm_selection_A(end), harm_selection_B(end), ...
    res_group.r_obs, res_group.p_perm, ...
    sum([res_ind.p_perm]<.05), length(subjects)));


if fig_save
    saveas(f, [fname_save_base, '_corrRdmsIndRho.png']);
end


%% compare rho with best model

r_A = [res_A.res.r]';
r_B = [res_B.res.r]';

z_A = r_to_z(r_A);
z_B = r_to_z(r_B);

c = [135, 175, 201]/255;

f = figure('pos', [1517 809 144 193], 'color', 'white');
pnl = panel(f);
ax = pnl.select();

s = scatter(repmat(0, 1, length(z_A)), z_A, 60, c, 'filled');
hold on ;
s = scatter(repmat(1, 1, length(z_B)), z_B, 60, c, 'filled');
plot(repmat([0;1],1,length(z_A)), [z_A, z_B]', '-', 'linew', 2, 'color', c);
ax.XLim = [-0.5, 1.5];

ax.XTick = [0, 1];
ax.XTickLabel = {sprintf('%s h-%d', measure, harm_selection_A(end)),
    sprintf('%s h-%d', measure, harm_selection_B(end))};
pnl.margin = [10,10,10,10];

[h, p, ci, stats] = ttest(z_A, z_B);

bf10 = bf.ttest(z_A, z_B);   % paired samples BF

txt = text(0.5, 1.1*max([z_A; z_B]), ...
    sprintf('r_A=%.2f r_B=%.2f\nBF10=%.2f, \np=%.2g', ...
    mean(r_A), mean(r_B), bf10, p));
txt.HorizontalAlignment = 'center';

disp(sprintf(['\n\nttest %s h-%d - h-%d (Fisher transformed), p = %.4f' ...
    '\nBF10 = %.4f'], ...
    measure, harm_selection_A(end), harm_selection_B(end), p, bf10));

if fig_save
    saveas(f, sprintf('%s_ttestRhoFisher.png', fname_save_base));
end


%% distribution of correlations across categorical model boundaries

r_A_distr = [res_A.res.r_all_models]';
r_B_distr = [res_B.res.r_all_models]';

z_A_distr = r_to_z(r_A_distr);
z_B_distr = r_to_z(r_B_distr);

r = diag(corr(z_A_distr', z_B_distr', 'type', 'pearson'));
z = r_to_z(r);

c_group_A = [31, 119, 180] / 255;      % Blue (for mean)
c_ind_A = [174, 199, 232] / 255;    % Light blue (for std)

c_group_B = [255, 127, 14] / 255;      % Orange (for mean)
c_ind_B = [255, 187, 120] / 255;    % Light orange (for std)

f = figure('pos', [476 675 358 191], 'color', 'white');
pnl = panel(f);
pnl.pack('h', [86, 14]);

ax = pnl(1).select();

tested_bounds = res_A.res(1).all_models_info.mod_id';

mu_A = mean(r_A_distr);
sem_A = std(r_A_distr, [], 1) ./ sqrt(size(r_A_distr, 1));

mu_B = mean(r_B_distr);
sem_B = std(r_B_distr, [], 1) ./ sqrt(size(r_B_distr, 1));

x = [tested_bounds, fliplr(tested_bounds)];
y = [mu_B + sem_B, fliplr(mu_B - sem_B)];
fill(x, y, c_ind_B, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
plot(tested_bounds, mu_B, 'LineWidth', 3.5, 'Color', c_group_B);
hold on;

x = [tested_bounds, fliplr(tested_bounds)];
y = [mu_A + sem_A, fliplr(mu_A - sem_A)];
fill(x, y, c_ind_A, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
plot(tested_bounds, mu_A, 'LineWidth', 3.5, 'Color', c_group_A);
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
ax.YLim = [-1, 4];
% ax.YTick = sort([min(ax.YTick), 0, 4]);
ylabel('Fisher-transformed r');

[h, p, ci, stats] = ttest(z);

txt = text(0, max(z)*1.2, sprintf('r = %.2f\np = %.2g', mean(r), p));
txt.HorizontalAlignment = 'center';

if fig_save
    saveas(f, [fname_save_base, '_distrCorr.png']);
end


disp(sprintf('\n\ncorrelation distribution rho %s h-%d - h-%d, r = %.4f p = %.10f', ...
    measure, harm_selection_A(end), harm_selection_B(end), mean(r), p));
