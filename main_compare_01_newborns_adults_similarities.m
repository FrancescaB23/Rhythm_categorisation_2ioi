% script to compare responses results from newborns and adults
% newborns: grand-average results
% adults: individual participants

clc
clear all
close all


%% parameters
par_adults = get_par('experiment', '2ioi_conga');
par_newborns = get_par('experiment', 'newborn');

ga_newborns = '19'; 

harm_selections = [1:6];
feat_snr = false;

roi_name = 'front';
avg_chans = true;

feat_resp = 'mX'; 

feat_stim = 'ratio';

n_perm_rsa_group = 10000;

fig_save = false;
run_stats = false;


%%
harm_str = vec_to_str(harm_selections, 'format','%d', 'sep','-'); 

if feat_snr
    snr_str = sprintf('%d-%d', par_adults.snr_bins(1), par_adults.snr_bins(2));
else
    snr_str = 'none';
end

save_path = fullfile(par_adults.deriv_path, 'comparisons', ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)));
mkdir(save_path);


%% load results
% load results adults
res_adults = load(fullfile( ...
    par_adults.summary_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf(['response-lf_feat-%s_harm-%s_snr-%s_' ...
    'roi-%s_avgChans-%s_cat-2-3-4-5_data.mat'], ...
    feat_resp, harm_str, snr_str, roi_name, jsonencode(avg_chans))));


% load results newborns
res_newborns = load(fullfile( ...
    par_newborns.results_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf('sub-grand-%s', ga_newborns), ...
    sprintf('sub-grand-%s_feat-%s_method-pearson_nCateg-2-3-4-5_perm.mat', ...
    ga_newborns, feat_resp)));

res_newborns.rdm = load(fullfile( ...
    par_newborns.results_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf('sub-grand-%s', ga_newborns), ...
    sprintf('sub-grand-%s_feat-%s_method-pearson_rdm.mat', ...
    ga_newborns, feat_resp)));

% rearrange results all models for ascending value of bond1, 2, 3, 4
[res_newborns.all_models_info, idx] = ...
    sortrows(res_newborns.all_models_info, ...
    {'n_categ', 'bound1', 'bound2', 'bound3', 'bound4'});
res_newborns.all_models_info.mod_id = [1 : length(res_newborns.all_models_info.n_categ)]';


%% sanity check
% check that the info about all models are ordered the same way across
% adults and newborns

vars_to_check = {'n_categ', 'bound1', 'bound2', 'bound3', 'bound4', 'bound5'};

if isequaln(res_adults.res(1).all_models_info{:, vars_to_check}, ...
        res_newborns.all_models_info{:, vars_to_check})
    disp('order models info match between groups');
else
    error('order models info do not match between groups');
end


%% save fname
fname_save_base = fullfile(save_path, ...
    sprintf('response-lf_feat-%s_featStim-%s_harm-%s_snr-%s_roi-%s_avgChans-%s_adultsVsNewborns', ...
    feat_resp, ...
    feat_stim, ...
    harm_str, ...
    snr_str, ...
    roi_name, ...
    jsonencode(avg_chans)));


%% stimuli info and rdm
load(fullfile(par_adults.deriv_path, 'stimuli', 'stim.mat'));
rdm_stim = rsa_stim_lf(par_newborns, [], 'feat', feat_stim);


%% correlate RSMs
% group-level permutation test
if exist([fname_save_base, '_permGroup.mat']) && ~run_stats
    clear res_group
    load([fname_save_base, '_permGroup.mat']);
else
    res_group = [];
    [res_group.r_obs, res_group.r_null_dist, res_group.p_perm] = ...
        perm_rdm_corr_group_empirical_model(...
        repmat({res_newborns.rdm.rdm}, 1, size({res_adults.res.rdm}, 2)),  ...
        {res_adults.res.rdm}, ...
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
    for i_sub=1:length(res_adults.res)
        i_sub
        [res_ind(i_sub).r_obs, res_ind(i_sub).p_perm, res_ind(i_sub).r_null_dist] = ...
            perm_rdm_corr(...
            res_newborns.rdm.rdm, ...
            res_adults.res(i_sub).rdm, ...
            n_perm_rsa_group, ...
            'parfor', true, ...
            'covariate', rdm_stim);
    end
    save([fname_save_base, '_permInd.mat'], 'res_ind');
end


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
       [fname_save_base, '_corrRdmsIndRho.png']); 
end


%% plot RSMs (individual adults)
f = figure('pos', [266 564 941 247], 'color', 'white'); % for 2 rows
pnl = panel(f);

plot_row = 2;
plot_col = ceil(length(res_adults.res)/plot_row);

pnl.pack('h', ones(1, plot_col)*floor(100/plot_col));

for i_h = [1:plot_col]
    pnl(i_h).pack('v', plot_row);
end

pnl.de.margin = 0;

count_row = 1;
count_col = 1;

for i_sub = 1:length(res_adults.res)

    ax = pnl(count_col, count_row).select();

    imagesc(ax, res_adults.res(i_sub).rdm);
%     title(sprintf('sub-%03d', res_adults.res(i_sub).sub));

    if ~mod(count_col, plot_col)
        count_col = 1;
        count_row = count_row + 1;

    else
         count_col = count_col + 1;
    end

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

pnl.de.margin = [0, 0, 0, 8];
pnl.margin = [20, 2, 2, 15];

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_RSMs_adults.png']);
end


%% RSM grand average newborn
f = figure('pos',[25 485 182 151], 'color', 'white');
imagesc(res_newborns.rdm.rdm);
title(sprintf('sub-grand-%s', ga_newborns));
axis off
axis image

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, sprintf('_RSM_newborn_GA-%s.png', ga_newborns)]);
end


%% stimulus RSM
f = figure('pos',[208 486 181 150], 'color', 'white');
imagesc(rdm_stim);
title(sprintf('stim-%s', feat_stim));
axis off
axis image

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, sprintf('_RSM_stim-%s.png', feat_stim)]);
end


%% distribution of correlation coefficients across categorical models
r_adults_distr = [res_adults.res.r_all_models]';
r_newborns_distr = res_newborns.all_models_info.r';

% Fisher's transformation
z_adults_distr = r_to_z(r_adults_distr);
z_newborns_distr = r_to_z(r_newborns_distr);

% correlate
for i_sub = 1:size(z_adults_distr, 1)

    r(i_sub) = corr(z_adults_distr(i_sub, :)', z_newborns_distr', 'type', 'pearson');
end

% Fisher's transformation
z = r_to_z(r);

% plot
z_mu = mean(z); 
z_sem = std(z) / sqrt(length(z)); 
z_ci = norminv(1 - 0.025) * z_sem; 

% f = figure('pos', [390 481 404 190], 'color', 'white'); % [1 90 435 220]
f = figure('pos', [476 675 358 191], 'color', 'white');

pnl = panel(f); 
pnl.pack('h', [86, 14]); 

ax = pnl(1).select();

tested_bounds = res_adults.res(1).all_models_info.mod_id';
mu_adults = mean(r_adults_distr);
sem_adults = std(r_adults_distr, [], 1) ./sqrt(size( ...
    r_adults_distr, 1));

plot(tested_bounds, mu_adults, 'LineWidth', 3);
hold on;
x = [tested_bounds, fliplr(tested_bounds)];
y = [mu_adults + sem_adults, fliplr(mu_adults - sem_adults)];
fill(x, y, 'b', 'FaceAlpha', 0.1, 'EdgeColor', 'none');
hold on;
plot(tested_bounds, r_newborns_distr, 'r', 'LineWidth', 3);
ax.XAxis.TickLength = [0,0]; 
ax.YAxis.TickLength = [0,0]; 
ax.YAxis.LineWidth = 2; 
ax.XAxis.LineWidth = 2; 
xlabel('categorical model');
ylabel('rho');
ylim([-.6 .65])


ax = pnl(2).select(); 
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
ax.YLim = [-1, 2]; 
ylabel('Fisher-transformed r');

[h, p, ci, stats] = ttest(z);     

txt = text(0, 1.8, sprintf('r = %.2f\np = %.2g', mean(r), p)); 
txt.HorizontalAlignment = 'center'; 

% if fig_save
%     print(f, '-dpng', '-painters', '-r600', ...
%         [fname_save_base, '_corrDistrRho.png']);
% end


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
ax.YLim = [min(z)-.25, max(z)+.25]; 
ax.YTick = [-.5, 0, 1];

ax.FontSize = 14; 
ax.FontName = 'Arial';
ylabel(ax, 'ρ (Fisher-transformed)', 'FontName', 'Arial', 'FontSize', 14);

% ylabel('Fisher-transformed r');

[h, p, ci, stats] = ttest(z);     

txt = title(ax, sprintf('r = %.2f\np = %.2g', mean(r), p)); 
txt.FontSize = 10;
txt.FontWeight = 'normal';

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_corrDistrRho.png']);
end

%% generate all possible categorical models (for plotting)
[models, all_models_info] = generate_categ_models_5(10, [2, 3, 4, 5], ...
    'n_mini_categ_skip', par_adults.n_mini_categ_skip);


%% plot correlation coefficients - descending order (by newborns)
% plot correlation coefficients across the 33 categorical models for eeg
% adults - newborns as dots, ordered by descending rho for newborn results

% sort r_newborns_distr by descending correlation coefficient
[r_newborns_distr_desc_nb, i_desc_nb] = sort(r_newborns_distr, 'descend');

r_adults_distr_desc_nb = zeros(size(r_adults_distr));

for i_sub = 1:size(r_adults_distr, 1)
    
    r_adults_distr_desc_nb(i_sub, :) = r_adults_distr(i_sub, i_desc_nb);

end


f = figure('Position', [1 496 629 370]); % 3 tiles [1 495 993 371] [1 495 1512 371] 
t = tiledlayout(3, length(models), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(t, [1 length(models)]);

for i_sub = 1:size(r_adults_distr_desc_nb, 1)
    h = scatter(ax1, ...
        (1:(length(models)))-.15, r_adults_distr_desc_nb(i_sub, :), ...
        70, [64 176 166]/255, 'filled');
    h.MarkerFaceAlpha = 0.60;
    hold on;

end

plot(ax1, (1:(length(models)))+0.15, r_newborns_distr_desc_nb, ...
    '^','MarkerFaceColor',[0 0 0], 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 12);
hold on; 
plot(ax1, (1:(length(models)))+0.15, r_newborns_distr_desc_nb, ...
    '^','MarkerFaceColor',[225 190 106]/255, 'MarkerEdgeColor', 'none', ...
    'MarkerSize', 10);

xlim(ax1, [0.5 (length(models))+0.5]);
title(ax1,sprintf('adults and newborns %s, descending by newborns', feat_resp));

labels = strcat('c', string(all_models_info.n_categ));
labels = labels(i_desc_nb);
xticks(ax1, 1:length(models));
xticklabels(ax1, labels);
% ylabel(ax1,'rho');
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
    imagesc(ax, models{i_desc_nb(i)});
    
    axis(ax,'off','square');
end

% index model
for i = 1:length(models)
    ax = nexttile(t);  % goes to row 3
    text(0.5, 0.5, sprintf('%d', i_desc_nb(i)), ...
         'HorizontalAlignment', 'center', ...
         'VerticalAlignment', 'middle', ...
         'FontSize', 10);
    axis(ax, 'off');
end

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_corrDistrDescNewborn.png']);
end


% plot models only
f = figure('Position', [1 771 1512 95]); % 3 tiles

n = length(models);
gap_units = 3;      % width of gap
tile_units = 10;    % width of real tile
cols = n*tile_units + (n-1)*gap_units;

t = tiledlayout(1, cols, ...
    'TileSpacing', 'none', 'Padding', 'compact');

col_start = 1;

for i = 1:n
    ax = nexttile(t, [1 tile_units]);   % real tile - span many columns
    imagesc(ax, models{i_desc_nb(i)});
    axis(ax,'off','square');
    
    col_start = col_start + tile_units;

    if i < n
        nexttile(t, [1 gap_units]);     % spacer tile -very thin gap
        axis off;
        col_start = col_start + gap_units;
    end
end


if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_corrDistrDescNewborn_models.png']);
end


%% Other figures
% plot correlation coefficients newborns and adults

f = figure('pos', [477 661 276 205], 'color', 'white'); % [827 506 256 190]
pnl = panel(f);
ax = pnl.select();

h = scatter(ax, 1.5, res_newborns.r_obs, 220, ... %140
    [0 0 0], '^', 'filled');
hold(ax, 'on');
h = scatter(ax, 1.5, res_newborns.r_obs, 150, ... %100
    [225, 190, 106] / 255, '^',  'filled');

hold(ax, 'on');

c_group = [38, 105, 99] / 255;
c_ind = [64, 176, 166] / 255;
r_all = [res_adults.res.r];
p_all = [res_adults.res.p_perm];

h = bar(ax, 4, mean(r_all), 'FaceColor', 'none', ...
        'EdgeColor', c_group, 'LineWidth', 3); 

hold(ax, 'on');

idx_signif = p_all < 0.05; 
    
h_signif = scatter(ax, 4 + 0.5 * (rand(1, sum(idx_signif))-0.5), r_all(idx_signif), ...
        110, c_ind, 'filled'); %50
h_signif.MarkerFaceAlpha = 0.7; % = 0.79;
h_signif.MarkerEdgeAlpha = 0; % = 0.79;

h_nonsignif = scatter(ax, 4 + 0.5 * (rand(1, sum(~idx_signif))-0.5), r_all(~idx_signif), ...
        110, c_ind); 
h_nonsignif.MarkerFaceAlpha = 0.79; % = 0.79;
h_nonsignif.MarkerEdgeAlpha = 0.79; % = 0.79;

sem = std(r_all) / sqrt(length(r_all)); 
ci = sem * norminv(1 - 0.025); 

plot(ax, [4,4], [mean(r_all)-ci, mean(r_all)+ci], ...
    'color', c_group, 'LineWidth', 4); 

box(ax, 'off'); 
ax.XLim = [1-0.7, 4.7]; 
ax.XTick = []; 
ax.YLim = [0, 1]; 
ax.YTick = [0, 1];
ax.FontSize = 16;
ax.FontName = 'Arial';
ax.YLabel.String = 'ρ'; 
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.TickLength = [0,0];  

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_corrRho_newborns_adults.png']);
end


%% plot overlay best theoretical model 
% adults
f = plot_rdms_overlay({res_adults.res.best_model}, 'alpha', 0.1);
title(sprintf('EEG ind adults %s', feat_resp));

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_overlay-best-model_adults.png']);
end

% grand-average newborns
f = plot_rdms_overlay({res_newborns.best_model}, 'alpha', 0.1);
title(sprintf('EEG GA newb %s', feat_resp));

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        [fname_save_base, '_overlay-best-model_newborns.png']);
end
