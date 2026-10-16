% script to compare responses results from eeg and tap

clear 
close all
clc


%% parameters
par = get_par('experiment', '2ioi_conga'); 
subjects = [num2cell([2:7 9:15 17:20])];

feat_name_resp = 'mX'; 
feat_name_stim = 'ratio';

feat_snr = false; 
roi_name = 'front';         
avg_chans = true; 
harm_selection = [1:6]; 

n_perm_rsa_group = 10000;

run_stats = false; 
fig_save = false;


%% colors for plotting
c_group_EEG = [38, 105, 99] / 255;
c_ind_EEG = [64, 176, 166] / 255;

c_group_TAP = [136, 34, 85] / 255;
c_ind_TAP = [204, 102, 119] / 255;


%%
if feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end

harm_str = vec_to_str(harm_selection, 'format','%d', 'sep','-'); 

save_path = fullfile(par.deriv_path, 'comparisons', 'eeg_vs_tap_mX');
mkdir(save_path); 


%% load data for EEG mX
response = 'lf';

subdirs = fullfile(sprintf('response-%s_stim-%s', response, feat_name_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)));

fname = fullfile(par.summary_path, subdirs, ...
    sprintf('response-%s_feat-%s_harm-%s_snr-%s_roi-%s_avgChans-%s_cat-2-3-4-5_data.mat', ...
    response, ...
    feat_name_resp, ...
    harm_str, ...
    snr_str, ...
    roi_name, ...
    jsonencode(avg_chans)));

res_EEG_mX = load(fname);


%% load data for TAP mX
response = 'tap-force';

subdirs = fullfile(sprintf('response-%s_stim-%s', response, feat_name_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str));

fname = fullfile(par.summary_path, subdirs, ...
    sprintf('response-%s_feat-%s_harm-%s_snr-%s_cat-2-3-4-5_data.mat', ...
    response, ...
    feat_name_resp, ...
    harm_str, ...
    snr_str));

res_TAP_mX = load(fname);


%% sanity check
if isequal([res_EEG_mX.res.sub], [res_TAP_mX.res.sub], [subjects{:}])
    disp('subject order across the two measures match, continuing the script');
else
    error('subject order across the two measures do not match, check!');
end


%% save fname 
fname_save_base = fullfile(save_path, ...
    sprintf(['response-eeg-versus-tap_feat-%s_featStim-%s_harm-%s_' ...
    'snr-%s_roi-%s_avgChans-%s'], ...
    feat_name_resp, ...
                     feat_name_stim, ...
                     harm_str, ...
                     snr_str, ...
                     roi_name, ...
                     jsonencode(avg_chans)));


%% correlate RDMs based on EEG and TAP
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
        {res_EEG_mX.res.rdm}, ...
        {res_TAP_mX.res.rdm},  ...
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
            res_EEG_mX.res(i_sub).rdm, ...
            res_TAP_mX.res(i_sub).rdm, ...
            n_perm_rsa_group, ...
            'parfor', true, ...
            'covariate', rdm_stim);
    end
    save([fname_save_base, '_permInd.mat'], 'res_ind');
end


disp(sprintf(['\n\npartial correlation, rho = %.4f, p = %.8f\n' ...
    'significant in %d participants out of %d (p<0.05)'], ...
    res_group.r_obs, res_group.p_perm, ...
    sum([res_ind.p_perm]<.05), length(subjects)));

%% plot partial correlation RSMs results
f = figure('pos', [827 508 153 189], 'color', 'white'); %[827 508 147 190]
pnl = panel(f);
ax = pnl.select();
plot_r_bar([res_ind.r_obs], [res_ind.p_perm],'y_maxlim', 1, 'ax', ax);
ax.YLim = [min([res_ind.r_obs])-.1, 1];
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.XAxisLocation = 'origin'; 
ax.YTick = [0, 1]; 
ax.FontSize = 14; 
ax.FontName = 'Arial';
ylabel(ax, 'ρ', 'FontName', 'Arial', 'FontSize', 14);
txt = text(0, 1, sprintf('r=%.2f\np=%.1g', res_group.r_obs, res_group.p_perm));
txt.HorizontalAlignment = 'center';

if fig_save
   print(f, '-dpng', '-painters', '-r600', ...
       [fname_save_base, '_corrRSMs.png']); 
end


%% plot individual RSM for EEG and TAP (individual)
% format long, with results
f = figure('pos', [1 84 575 782], 'color', 'white');
pnl = panel(f);

pnl.pack('v', [ones(1, length(res_EEG_mX.res))*floor(90 / length(res_EEG_mX.res))]);

for i_h = 1:size(subjects, 2)
    pnl(i_h).pack('h', 5);
end

% Adjust column margins to add space between columns 1-2 and 3-4
pnl.de.marginleft = 5;
pnl.de.marginright = 5;

for i_h = 1:size(subjects, 2)

    % subject code
    ax = pnl(i_h, 1).select();
    axis(ax, 'off');
    text(ax, 0.5, 0.5, sprintf('sub-%03d', res_EEG_mX.res(i_h).sub), ...
        'HorizontalAlignment', 'center', ...
        'VerticalAlignment', 'middle', ...
        'FontName', 'Arial', ...
        'FontSize', 12);

    % EEG RSM
    ax = pnl(i_h, 2).select();
    imagesc(ax, res_EEG_mX.res(i_h).rdm);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');
    ax.XTick = [];
    ax.YTick = [];
    % make sure the y-axis direction is reverse (the default differes depending
    % on whether the axes are obtained from the builtin subplot function or
    % from panel)
    ax.YDir = 'reverse';
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';


    % EEG best 2-cat model RSM
    ax = pnl(i_h, 3).select();
    imagesc(ax, res_EEG_mX.res(i_h).best_model);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');
    ax.XTick = [];
    ax.YTick = [];
    ax.YDir = 'reverse';
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';
    

    % TAP RSM
    ax = pnl(i_h, 4).select();
    imagesc(ax, res_TAP_mX.res(i_h).rdm);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');
    ax.XTick = [];
    ax.YTick = [];
    ax.YDir = 'reverse';
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';

    % TAP best 2-cat model RSM
    ax = pnl(i_h, 5).select();
    imagesc(ax, res_TAP_mX.res(i_h).best_model);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');
    ax.XTick = [];
    ax.YTick = [];
    ax.YDir = 'reverse';
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';

end

pnl.de.margin = [5, 5, 5, 5];
pnl.margin = [20, 2, 2, 15];
for i_h = 1:size(subjects, 2)
    pnl(i_h, 1).marginright = 30;  % Add space after column 1
    pnl(i_h, 3).marginright = 30;  % Add space after column 3
end

% Add centered titles spanning columns 2-3 and 4-5 for each row
for i_h = 1:size(subjects, 2)
    % Get positions of columns 2 and 3 for EEG title
    pos2 = pnl(i_h, 2).position;
    pos3 = pnl(i_h, 3).position;
    center_x_EEG = (pos2(1) + pos3(1) + pos3(3)) / 2;
    title_y = pos2(2) + pos2(4) + 0.00001;  % Slightly above the plots
    
    % Add EEG title centered over columns 2-3
    annotation('textbox', [center_x_EEG - 0.15, title_y, 0.3, 0.02], ...
        'String', sprintf('r=%.4f p_B=%.4f', res_EEG_mX.res(i_h).r, res_EEG_mX.res(i_h).p_perm), ...
        'HorizontalAlignment', 'center', ...
        'EdgeColor', 'none', ...
        'FontName', 'Arial', ...
        'FontWeight', 'bold', ...
        'FontSize', 10);
    
    % Get positions of columns 4 and 5 for TAP title
    pos4 = pnl(i_h, 4).position;
    pos5 = pnl(i_h, 5).position;
    center_x_TAP = (pos4(1) + pos5(1) + pos5(3)) / 2;
    
    % Add TAP title centered over columns 4-5
    annotation('textbox', [center_x_TAP - 0.15, title_y, 0.3, 0.02], ...
        'String', sprintf('r=%.4f p_B=%.4f', ...
        res_TAP_mX.res(i_h).r, res_TAP_mX.res(i_h).p_perm), ...
        'HorizontalAlignment', 'center', ...
        'EdgeColor', 'none', ...
        'FontName', 'Arial', ...
        'FontWeight', 'bold', ...
        'FontSize', 10);
end

if fig_save
   print(f, '-dpng', '-painters', '-r600', ...
       [fname_save_base, '_individualResults.png']); 
end


%% distribution of correlations across categorical model boundaries
r_EEG_mX_distr = [res_EEG_mX.res.r_all_models]';
r_TAP_mX_distr = [res_TAP_mX.res.r_all_models]';

z_EEG_mX_distr = r_to_z(r_EEG_mX_distr);
z_TAP_mX_distr = r_to_z(r_TAP_mX_distr);

r = diag(corr(z_EEG_mX_distr', z_TAP_mX_distr', 'type', 'pearson'));
z = r_to_z(r);

f = figure('pos', [476 675 358 191], 'color', 'white');
pnl = panel(f);
pnl.pack('h', [86, 14]); 

ax = pnl(1).select();
tested_bounds = res_EEG_mX.res(1).all_models_info.mod_id';

mu_EEG = mean(r_EEG_mX_distr);
sem_EEG = std(r_EEG_mX_distr, [], 1) ./ sqrt(size(r_EEG_mX_distr, 1));

mu_TAP = mean(r_TAP_mX_distr);
sem_TAP = std(r_TAP_mX_distr, [], 1) ./ sqrt(size(r_TAP_mX_distr, 1));

x = [tested_bounds, fliplr(tested_bounds)];
y = [mu_EEG + sem_EEG, fliplr(mu_EEG - sem_EEG)];
fill(x, y, c_ind_EEG, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
plot(tested_bounds, mu_EEG, 'LineWidth', 3.5, 'Color', c_group_EEG);
hold on;

x = [tested_bounds, fliplr(tested_bounds)];
y = [mu_TAP + sem_TAP, fliplr(mu_TAP - sem_TAP)];
fill(x, y, c_ind_TAP, 'FaceAlpha', 0.5, 'EdgeColor', 'none');
hold on;
plot(tested_bounds, mu_TAP, 'LineWidth', 3.5, 'Color', c_group_TAP);
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
% ax.YLim = [-1, 3]; 
ylabel('Fisher-transformed r');

[h, p, ci, stats] = ttest(z);

txt = text(0, 3, sprintf('r = %.2f\np = %.2g', mean(r), p));
txt.HorizontalAlignment = 'center';


%% figure correlation distribution correlation coefficients
f = figure('pos', [827 508 120 189], 'color', 'white'); %[827 508 147 190] [827 508 153 189]
pnl = panel(f);
ax = pnl.select();
pnl.margin = [18 5 2 12];
c = [135, 175, 201] / 255; 
plot(0 + 0.2*rand(1,length(z))-0.1, z, 'o', 'color', c, 'markerfacecolor', c, 'markersize', 5); 
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
ax.YLim = [min(z)-.25, max(z)+.35]; 
ax.YTick = [-.5, 0, 1, 2];

ax.FontSize = 14; 
ax.FontName = 'Arial';
ylabel(ax, 'ρ (Fisher-transformed)', 'FontName', 'Arial', 'FontSize', 14);

[h, p, ci, stats] = ttest(z);     

txt = title(ax, sprintf('r = %.2f\np = %.2g', mean(r), p)); 
txt.FontSize = 10;
txt.FontWeight = 'normal';

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_corrDistrRho.png']);
end


%% plot correlation coefficients - EEG TAP - descending order by adults
% by TAP or by EEG (average adults)

order_by = 'tap'; % 'eeg' 'tap'

switch order_by
    case 'eeg'
        [~, i_desc_ad] = sort(mean(r_EEG_mX_distr), 'descend');
    case 'tap'
        [~, i_desc_ad] = sort(mean(r_TAP_mX_distr), 'descend');
end

% reorder rho for EEG and TAP
rho_cat_EEG_desc_ad = zeros(size(r_EEG_mX_distr));
rho_cat_TAP_desc_ad = zeros(size(r_TAP_mX_distr));

for i_sub = 1:size(r_EEG_mX_distr, 1)
    rho_cat_EEG_desc_ad(i_sub, :) = r_EEG_mX_distr(i_sub, i_desc_ad);
    rho_cat_TAP_desc_ad(i_sub, :) = r_TAP_mX_distr(i_sub, i_desc_ad);
end

% re-generate models (for plotting)
n_categ = [2, 3, 4, 5];

[models, all_models_info] = generate_categ_models_5(par.n_cond, n_categ, ...
    'n_mini_categ_skip', par.n_mini_categ_skip);

% sanity check: regenerated models order match the order of the results
cols_check = {'n_categ', 'bound1', 'bound2', 'bound3', 'bound4'};

if isequaln(res_EEG_mX.res(1).all_models_info{:, cols_check}, all_models_info{:, cols_check}) && ...
        isequaln(res_TAP_mX.res(1).all_models_info{:, cols_check}, all_models_info{:, cols_check})
    disp('model definitions match across eeg, tap, and re-generated models, continuing the script');
else
    error('model definitions do not match the regenerated ones, check!');
end


f = figure('Position', [1 496 629 370]); % 3 tiles [1 495 993 371] [1 495 1512 371] 
t = tiledlayout(3, length(models), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(t, [1 length(models)]);

for i_sub = 1:size(rho_cat_EEG_desc_ad, 1)
    h = scatter(ax1, ...
        (1:(length(models)))-.15, rho_cat_EEG_desc_ad(i_sub, :), ...
        70, c_ind_EEG, 'filled');
    h.MarkerFaceAlpha = 0.60;
    hold on;

end

for i_sub = 1:size(rho_cat_TAP_desc_ad, 1)

    h = scatter(ax1, ...
        (1:(length(models)))+.15, rho_cat_TAP_desc_ad(i_sub, :), ...
        70, c_ind_TAP, 'filled');
    h.MarkerFaceAlpha = 0.60;
    hold on;

end

xlim(ax1, [0.5 (length(models))+0.5]);
title(ax1,sprintf('adults EEG and TAP %s, descending by rho adults %s', ...
    feat_name_resp, order_by));

labels = strcat('c', string(all_models_info.n_categ));
labels = labels(i_desc_ad);
xticks(ax1, 1:length(models));
xticklabels(ax1, labels);
ax1.TickLength = [0, 0]; 
ax1.YLim = [-1, 1];
ax1.YTick = [-1, 0, 1];
ax1.YAxis.LineWidth = 2; 
ax1.XAxis.LineWidth = 2; 
ax1.FontSize = 12; 
ax1.FontName = 'Arial';
ylabel(ax1, 'ρ', 'FontName', 'Arial', 'FontSize', 12);

for i = 1:length(models)
    ax = nexttile(t); % goes to row 2
    imagesc(ax, models{i_desc_ad(i)});
    
    axis(ax,'off','square');
end

% index model
for i = 1:length(models)
    ax = nexttile(t);  % goes to row 3
    text(0.5, 0.5, sprintf('%d', i_desc_ad(i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'FontSize', 10);
    axis(ax, 'off');
end

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, sprintf('_corrDistrDescBy-%s.png', order_by)]);
end


% plot models only - ordered 
f = figure('Position', [1 771 1512 95]); % 3 tiles
n = length(models);
gap_units = 3;      % width of gap
tile_units = 10;    % width of real tile
cols = n*tile_units + (n-1)*gap_units;
t = tiledlayout(1, cols, ...
    'TileSpacing', 'none', 'Padding', 'compact');

col_start = 1;

for i = 1:n
    ax = nexttile(t, [1 tile_units]);   % span many columns
    imagesc(ax, models{i_desc_ad(i)});
    axis(ax,'off','square');
    
    col_start = col_start + tile_units;

    if i < n
        nexttile(t, [1 gap_units]);     % very thin gap
        axis off;
        col_start = col_start + gap_units;
    end
end

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, sprintf('_corrDistrDescBy-%s_models.png', order_by)]);
end


%% Other figures
% plot correlation coefficients eeg and tap
f = figure('pos', [477 661 276 205], 'color', 'white'); % [827 506 256 190]
pnl = panel(f);
ax = pnl.select();
hold(ax, 'on');

data_plot = {res_EEG_mX, res_TAP_mX};
x_pos  = [1.5, 4];
labels = {'EEG', 'TAP'};
c_group = {c_group_EEG, c_group_TAP};
c_ind = {c_ind_EEG,   c_ind_TAP};

for i = 1:numel(data_plot)
    r_all = [data_plot{i}.res.r];
    p_all = [data_plot{i}.res.p_perm];

    bar(ax, x_pos(i), mean(r_all), 'FaceColor', 'none', ...
        'EdgeColor', c_group{i}, 'LineWidth', 3, 'BarWidth', 0.9);

    idx_signif = p_all < 0.05;

    h_signif = scatter(ax, x_pos(i) + 0.5 * (rand(1, sum(idx_signif)) - 0.5), ...
        r_all(idx_signif), 110, c_ind{i}, 'filled');
    h_signif.MarkerFaceAlpha = 0.7;
    h_signif.MarkerEdgeAlpha = 0;

    h_nonsignif = scatter(ax, x_pos(i) + 0.5 * (rand(1, sum(~idx_signif)) - 0.5), ...
        r_all(~idx_signif), 110, c_ind{i});
    h_nonsignif.MarkerEdgeAlpha = 0.79;

    sem = std(r_all) / sqrt(numel(r_all));
    ci  = sem * norminv(1 - 0.025);
    plot(ax, [x_pos(i), x_pos(i)], [mean(r_all) - ci, mean(r_all) + ci], ...
        'Color', c_group{i}, 'LineWidth', 4);
end

box(ax, 'off');
ax.XLim = [x_pos(1) - 1.2, x_pos(2) + 1.2];
ax.XTick = [];
% ax.YLim = [0, 1];
ax.YTick = [0, 1];
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.YLabel.String = 'ρ';
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.TickLength = [0, 0];

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_corrRho_eeg_tap.png']);
end


%% Individual RSMs
% plot individual RSM for stim, EEG and TAP (individual)
f = figure('pos', [1 593 1422 273], 'color', 'white');
pnl = panel(f);

pnl.pack('h', [12 ones(1, length(res_EEG_mX.res))*floor(90 / length(res_EEG_mX.res))]);

for i_h = 2:(size(subjects, 2)+1)
    pnl(i_h).pack('v', 2);
end

ax = pnl(1).select();
imagesc(ax, rdm_stim);

set(ax, 'YDir', 'normal');
axis(ax, 'square');
ax.XTick = [];
ax.YTick = [];
% make sure the y-axis direction is reverse (the default differes depending
% on whether the axes are obtained from the builtin subplot function or
% from panel)
ax.YDir = 'reverse';
ax.XAxis.Visible = 'off';
ax.YAxis.Visible = 'off';

for i_h = 2:(size(subjects, 2)+1)

    % EEG RSM
    ax = pnl(i_h, 1).select();
    imagesc(ax, res_EEG_mX.res(i_h-1).rdm);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');
    ax.XTick = [];
    ax.YTick = [];
    % make sure the y-axis direction is reverse (the default differes depending
    % on whether the axes are obtained from the builtin subplot function or
    % from panel)
    ax.YDir = 'reverse';
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';

    % TAP RSM
    ax = pnl(i_h, 2).select();
    imagesc(ax, res_TAP_mX.res(i_h-1).rdm);
    set(ax, 'YDir', 'normal');
    axis(ax, 'square');
    ax.XTick = [];
    ax.YTick = [];
    % make sure the y-axis direction is reverse (the default differes depending
    % on whether the axes are obtained from the builtin subplot function or
    % from panel)
    ax.YDir = 'reverse';
    ax.XAxis.Visible = 'off';
    ax.YAxis.Visible = 'off';

end

pnl.de.margin = [5, 5, 5, 5];
pnl.margin = [20, 2, 2, 15];

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_RSM_individual.png']);
end


%% plot overlay best theoretical model 
% eeg
f = plot_rdms_overlay({res_EEG_mX.res.best_model}, 'alpha', 0.1);
title(sprintf('EEG ind adults %s', feat_name_resp));

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_overlay-best-model_eeg.png']);
end

% tap
f = plot_rdms_overlay({res_TAP_mX.res.best_model}, 'alpha', 0.1);
title(sprintf('TAP ind adults %s', feat_name_resp));

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_overlay-best-model_tap.png']);
end