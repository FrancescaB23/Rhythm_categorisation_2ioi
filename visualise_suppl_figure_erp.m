%% Supplementary figure - ERP results
% script to create the supplementary figure with the ERP results
% showing:
% stimuli (one repetition pattern) overlaid (envelope)
% newborn: grand-average ERP, RSM, plot rho, best categorical model
% adults, eeg and tap: one individual participant ERP and RSM, group rho,
% overlay best categorical model

clc
clear
close all


%% parameters
par_adults = get_par('experiment', '2ioi_conga');
par_newborns = get_par('experiment', 'newborn');

ga_newborns = '19';
sample_sub_adult = 9;

feat_stim = 'ratio';

roi_name = 'front';
avg_chans = true;

lpf_cutoff = 10;

fig_save = true;


%%
fig_folder = fullfile(par_adults.deriv_path, 'figures', 'supplementary');


%% load stimuli
load(fullfile(par_adults.deriv_path, 'stimuli', 'stim.mat'));


%% plot stimuli - erp like
% stimuli wave file
for i_cond = 1:size(stim.s, 1)
    s(i_cond, :) = stim.s{i_cond}(1:stim.fs*stim.cycle_dur);
end

% stimuli envelope
stim_env = cellfun(@(x) abs(hilbert(x)), stim.s, 'uni', 0);
for i_cond = 1:size(stim.s, 1)
    stim_env_cycle(i_cond, :) = stim_env{i_cond}(1:stim.fs*stim.cycle_dur);
end

f = plot_erp_overlay_F(stim_env_cycle, stim.fs, 'ratios', stim.ratios);
ax = findobj(f, 'Type', 'axes');
ax.YTick = [0, 1];
ax.FontSize = 14;
ax.FontName = 'Arial';
ylabel(ax, 'Amplitude [a.u.]', 'FontName', 'Arial', 'FontSize', 15);

if fig_save
    fname_save = 'ERP_stim-env';
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_folder, fname_save), '.png']);
end


%% newborn results - time domain data for visualisation

% load data
load_path = fullfile(par_newborns.deriv_path, 'erp', ...
    sprintf('sub-grand-%s', ga_newborns));
fname = sprintf('sub-grand-%s_response-erp_lpf-%.0fHz_rr-m', ...
    ga_newborns, lpf_cutoff);
[header, data] = CLW_load(fullfile(load_path, fname));

% select (and average) channels
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
end
[header, data] = RLW_arrange_channels(header, data, roi_chans);

if avg_chans
    [header, data] = RLW_pool_channels(header, data, roi_chans, ...
        'mixed_channel_label', roi_name, ...
        'keep_original_channels', false);
end

% plot ERPs
f = plot_erp_overlay_F(squeeze(data), 1/header.xstep, ...
    'ratios', stim.ratios);
ax = findobj(f, 'Type', 'axes');
ax.YTick = [-.5, 0, .5];
ax.FontSize = 14;
ax.FontName = 'Arial';
ylabel(ax, 'Amplitude [μV]', 'FontName', 'Arial', 'FontSize', 15);

if fig_save
    fname_save = sprintf('ERP_eeg-newborn-erp_sub-grand-%s_roi-%s_avgChans-%s', ...
        ga_newborns, roi_name, jsonencode(avg_chans));
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_folder, fname_save), '.png']);
end


%% newborn results - RSA
% load results
subdirs = fullfile(par_newborns.results_path, ...
    sprintf('response-erp_stim-%s', feat_stim), ...
    sprintf('lfp-%dHz', lpf_cutoff), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf('sub-grand-%s', ga_newborns));
fname = sprintf('sub-grand-%s_feat-erp_method-pearson_', ga_newborns);

res_newborns = load(fullfile(subdirs, [fname, 'nCateg-2-3-4-5_perm.mat']))
res_newborns.rdm = load(fullfile(subdirs, [fname, 'rdm.mat']));

% plot RSM
f = figure('pos',[25 485 182 151], 'color', 'white');
imagesc(res_newborns.rdm.rdm);
title(sprintf('sub-grand-%s ERP', ga_newborns));
axis off
axis image

if fig_save
    fname_save = sprintf('ERP_eeg-newborn-rdm_sub-grand-%s_roi-%s_avgChans-%s', ...
        ga_newborns, roi_name, jsonencode(avg_chans));
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


% plot correlation coefficient
f = figure('pos', [1 698 139 172], 'color', 'white'); % [827 506 256 190] [477 661 276 205]
pnl = panel(f);
ax = pnl.select();

h = scatter(ax, 1, res_newborns.r_obs, 220, [0 0 0], '^', 'filled');
hold(ax, 'on');
h = scatter(ax, 1, res_newborns.r_obs, 150, [225, 190, 106] / 255, '^',  'filled');

box(ax, 'off');
ax.XTick = [];
ax.YLim = [0, 1];
ax.YTick = [0, 1];
ax.FontSize = 14;
ax.FontName = 'Arial';
% ax.YLabel.String = 'ρ';
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.TickLength = [0,0];
ylabel(ax, 'ρ', 'FontName', 'Arial', 'FontSize', 15);

if fig_save
    fname_save = sprintf('ERP_eeg-newborn-rho_sub-grand-%s_roi-%s_avgChans-%s', ...
        ga_newborns, roi_name, jsonencode(avg_chans));
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


% plot overlay best theoretical model - GA newborns
f = plot_rdms_overlay({res_newborns.best_model}, 'alpha', 0.1);
title('EEG GA newb ERP');

if fig_save
    fname_save = sprintf('ERP_eeg-newborn-best-model_sub-grand-%s_roi-%s_avgChans-%s', ...
        ga_newborns, roi_name, jsonencode(avg_chans));
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


%% adult EEG results - time domain for visualisation
% load data
sub_str = sub_num2str(sample_sub_adult);
load_path = fullfile(par_adults.deriv_path, 'erp',sub_str);
fname = sprintf('%s_response-erp_lpf-%.0fHz_rr-m', sub_str, lpf_cutoff);
[header, data] = CLW_load(fullfile(load_path, fname));

% select (and average) channels
if strcmp(roi_name, 'front')
    roi_chans = {
        'F1'
        'Fz'
        'F2'
        'FC1'
        'FCz'
        'FC2'
        'C1'
        'Cz'
        'C2'
        };
end
[header, data] = RLW_arrange_channels(header, data, roi_chans);

if avg_chans
    [header, data] = RLW_pool_channels(header, data, roi_chans, ...
        'mixed_channel_label', roi_name, ...
        'keep_original_channels', false);
end

% plot ERPs
f = plot_erp_overlay_F(squeeze(data), 1/header.xstep, 'ratios', stim.ratios);
ax = findobj(f, 'Type', 'axes');
ax.YTick = [-2, 0, 4];
ax.FontSize = 14;
ax.FontName = 'Arial';
ylabel(ax, 'Amplitude [μV]', 'FontName', 'Arial', 'FontSize', 15);

if fig_save
    fname_save = sprintf('ERP_eeg-adult-erp_%s_roi-%s_avgChans-%s', ...
        sub_str, roi_name, jsonencode(avg_chans));
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


%% adult EEG results - RSA
% load results
res_adults_eeg = load(fullfile(par_adults.summary_path, ...
    sprintf('response-erp_stim-%s', feat_stim), ...
    sprintf('lpf-%dHz', lpf_cutoff), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf('response-erp_lpf-%dHz_roi-%s_avgChans-%s_cat-2-3-4-5_data.mat', ...
    lpf_cutoff, roi_name, jsonencode(avg_chans))));


% plot RSM
idx_sample_sub_adult = find([res_adults_eeg.res.sub] == sample_sub_adult);
f = figure('pos',[25 485 182 151], 'color', 'white');
imagesc(res_adults_eeg.res(idx_sample_sub_adult).rdm);
title(sprintf('%s EEG ERP', sub_str));
axis off
axis image

if fig_save
    fname_save = sprintf('ERP_eeg-adult-rdm_%s_roi-%s_avgChans-%s', ...
        sub_str, roi_name, jsonencode(avg_chans));
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


% plot correlation coefficient
f = figure('pos', [1 698 139 172], 'color', 'white'); % [827 506 256 190] [477 661 276 205]
pnl = panel(f);
ax = pnl.select();

hold(ax, 'on');

c_group = [38, 105, 99] / 255;
c_ind = [64, 176, 166] / 255;
r_all = [res_adults_eeg.res.r];
p_all = [res_adults_eeg.res.p_perm];
h = bar(ax, 1, mean(r_all), 'FaceColor', 'none', ...
    'EdgeColor', c_group, 'LineWidth', 3);
hold(ax, 'on');
idx_signif = p_all < 0.05;
h_signif = scatter(ax, 1 + 0.5 * (rand(1, sum(idx_signif))-0.5), r_all(idx_signif), ...
    110, c_ind, 'filled'); %50
h_signif.MarkerFaceAlpha = 0.7; % = 0.79;
h_signif.MarkerEdgeAlpha = 0; % = 0.79;
h_nonsignif = scatter(ax, 1 + 0.5 * (rand(1, sum(~idx_signif))-0.5), r_all(~idx_signif), ...
    110, c_ind);
h_nonsignif.MarkerFaceAlpha = 0.79; % = 0.79;
h_nonsignif.MarkerEdgeAlpha = 0.79; % = 0.79;
sem = std(r_all) / sqrt(length(r_all));
ci = sem * norminv(1 - 0.025);
plot(ax, [1, 1], [mean(r_all)-ci, mean(r_all)+ci], ...
    'color', c_group, 'LineWidth', 4);

box(ax, 'off');
ax.YLim = [0, 1];
ax.YTick = [0, 1];
ax.XTick = [];
ax.FontSize = 14;
ax.FontName = 'Arial';
% ax.YLabel.String = 'ρ';
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.TickLength = [0,0];
ylabel(ax, 'ρ', 'FontName', 'Arial', 'FontSize', 15);

if fig_save
    fname_save = sprintf('ERP_eeg-adult-rho_group_roi-%s_avgChans-%s', ...
        roi_name, jsonencode(avg_chans));
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


% plot overlay best theoretical model - adults
f = plot_rdms_overlay({res_adults_eeg.res.best_model}, 'alpha', 0.1);
title('EEG ind adults ERP');

if fig_save
    fname_save = sprintf('ERP_eeg-adult-best-model_group_roi-%s_avgChans-%s', ...
        roi_name, jsonencode(avg_chans));
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


%% adult TAP results - time domain for visualisation
% load data
sub_str = sub_num2str(sample_sub_adult);
load_path = fullfile(par_adults.deriv_path,'erp', sub_str);
fname = sprintf('%s_response-tap-force-erp_lpf-%.0fHz', ...
    sub_str, lpf_cutoff);
[header, data] = CLW_load(fullfile(load_path, fname));

% plot ERPs
f = plot_erp_overlay_F(squeeze(mean(data, 2)), 1/header.xstep, ...
    'ratios', stim.ratios);
ax = findobj(f, 'Type', 'axes');
ax.YTick = [0];
ax.FontSize = 14;
ax.FontName = 'Arial';
ylabel(ax, 'Tap force [a.u.]', 'FontName', 'Arial', 'FontSize', 15);

if fig_save
    fname_save = sprintf('ERP_tap-adult-erp_%s', sub_str);
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


%% adult TAP results - RSA
% load results
res_adults_tap = load(fullfile(par_adults.summary_path, ...
    sprintf('response-tap-force-erp_stim-%s', feat_stim), ...
    sprintf('lpf-%dHz', lpf_cutoff), ...
    sprintf('response-tap-force-erp_lpf-%dHz_cat-2-3-4-5_data.mat', lpf_cutoff)));


% plot RSM
idx_sample_sub_adult = find([res_adults_tap.res.sub] == sample_sub_adult);
f = figure('pos',[25 485 182 151], 'color', 'white');
imagesc(res_adults_tap.res(idx_sample_sub_adult).rdm);
title(sprintf('%s TAP ERP', sub_str));
axis off
axis image

if fig_save
    fname_save = sprintf('ERP_tap-adult-rdm_%s', sub_str);
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


% plot correlation coefficient
f = figure('pos', [1 698 139 172], 'color', 'white'); % [827 506 256 190] [477 661 276 205]  % [2 665 169 205]
pnl = panel(f);
ax = pnl.select();

hold(ax, 'on');

c_group = [136, 34, 85] / 255;
c_ind = [204, 102, 119] / 255;
r_all = [res_adults_tap.res.r];
p_all = [res_adults_tap.res.p_perm];

h = bar(ax, 1, mean(r_all), 'FaceColor', 'none', ...
    'EdgeColor', c_group, 'LineWidth', 3);
hold(ax, 'on');
idx_signif = p_all < 0.05;
h_signif = scatter(ax, 1 + 0.5 * (rand(1, sum(idx_signif))-0.5), r_all(idx_signif), ...
    110, c_ind, 'filled'); %50
h_signif.MarkerFaceAlpha = 0.7; % = 0.79;
h_signif.MarkerEdgeAlpha = 0; % = 0.79;
h_nonsignif = scatter(ax, 1 + 0.5 * (rand(1, sum(~idx_signif))-0.5), r_all(~idx_signif), ...
    110, c_ind);
h_nonsignif.MarkerFaceAlpha = 0.79; % = 0.79;
h_nonsignif.MarkerEdgeAlpha = 0.79; % = 0.79;
sem = std(r_all) / sqrt(length(r_all));
ci = sem * norminv(1 - 0.025);
plot(ax, [1, 1], [mean(r_all)-ci, mean(r_all)+ci], ...
    'color', c_group, 'LineWidth', 4);
box(ax, 'off');
ax.XTick = [];
ax.YLim = [0, 1];
ax.YTick = [0, 1];
ax.FontSize = 14;
ax.FontName = 'Arial';
% ax.YLabel.String = 'ρ';
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.TickLength = [0,0];
ylabel(ax, 'ρ', 'FontName', 'Arial', 'FontSize', 15);


if fig_save
    fname_save = 'ERP_tap-adult-rho_group';
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end


% plot overlay best theoretical model - adults tap
f = plot_rdms_overlay({res_adults_tap.res.best_model}, 'alpha', 0.1);
title('TAP ind adults ERP');

if fig_save
    fname_save = 'ERP_tap-adult-best-model_group';
    print(f, '-dpng', '-painters', '-r600', [fullfile(fig_folder, fname_save), '.png']);
end