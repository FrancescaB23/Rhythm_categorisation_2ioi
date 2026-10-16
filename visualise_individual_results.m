% script to plot individual results for the figure describing the method

clc
clear all
close all


%% parameters
par = get_par('experiment', '2ioi_conga');

harm_selections = [1:6];
feat_snr = false;

roi_name = 'front';
avg_chans = true;

feat_resp = 'mX'; 

feat_stim = 'ratio';

fig_save = true;


%%
harm_str = vec_to_str(harm_selections, 'format','%d', 'sep','-'); 

if feat_snr
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2));
else
    snr_str = 'none';
end

figure_path = fullfile(par.deriv_path, 'figures');

if ~exist(figure_path)
    mkdir(figure_path);
end


%% load results
res_adults = load(fullfile( ...
    par.summary_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf(['response-lf_feat-%s_harm-%s_snr-%s_' ...
    'roi-%s_avgChans-%s_cat-2-3-4-5_data.mat'], ...
    feat_resp, harm_str, snr_str, roi_name, jsonencode(avg_chans))));


%% individual plots for visualisation
% distribution of correlation coefficients and RSM of one example
% participants for Figure 1

i_sub = 7; % 5 or 7

% isequal(res_adults.res(i_sub).r_all_models', r_adults_distr(i_sub, :))

% plot correlation coefficients
f = figure('Position', [1 745 279 118]); % [1 678 467 156]
t = tiledlayout(1, length(res_adults.res(i_sub).r_all_models), ...
    'TileSpacing', 'compact', 'Padding', 'compact');
ax1 = nexttile(t, [1 length(res_adults.res(i_sub).r_all_models)]);

h = scatter(ax1, ...
    (1:(length(res_adults.res(i_sub).r_all_models))), ...
    res_adults.res(i_sub).r_all_models, ...
    55, [64 176 166]/255, 'filled');
h.MarkerFaceAlpha = 1;

xlim(ax1, [0.5 (length(res_adults.res(i_sub).r_all_models))+0.5]);
title(ax1,sprintf('adults sub-%03d EEG %s', ...
    res_adults.res(i_sub).sub, feat_resp));

labels = [1:length(res_adults.res(i_sub).r_all_models)];
xticks(ax1, 1:length(res_adults.res(i_sub).r_all_models));
xticklabels(ax1, labels);
ax1.TickLength = [0, 0]; 
ax1.YLim = [-1, 1];
yticklabels(ax1, {});
ax1.YAxis.LineWidth = 2; 
ax1.XAxis.LineWidth = 2; 

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        fullfile(figure_path, ...
        sprintf('visualisation_sub-%03d_EEG_distr-rho.png', ...
        res_adults.res(i_sub).sub)));
end

%% plot RSM
f = figure('pos',[25 485 182 151], 'color', 'white');
imagesc(res_adults.res(i_sub).rdm);
title(sprintf('adults sub-%03d EEG %s', ...
    res_adults.res(i_sub).sub, feat_resp));
axis off
axis image

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        fullfile(figure_path, ...
        sprintf('visualisation_sub-%03d_EEG_RSM.png', ...
        res_adults.res(i_sub).sub)));
end

