%% script to create the figures for tap onset extraction and ITI ratios

clc
clear
close all


%% parameters
par = get_par('experiment', '2ioi_conga');

sub_tap_onset = 10;
cond_tap_onset = 10;
n_cycle_plot = 10;

response = 'tap-onsets';
feat_name_stim = 'ratio';

subjects = [2:7 9:15 17:20];

fig_save = true;


%%
fig_folder = fullfile(par.deriv_path, 'figures', 'supplementary');


%% plot tap onset extraction
% load stimuli
load(fullfile(par.deriv_path, 'stimuli', 'stim.mat')); 


% load tapping data
sub_str = sub_num2str(sub_tap_onset);
load_path = fullfile(par.deriv_path, 'tapping', sub_str);
fname = sprintf('%s_task-TAP_response-force', sub_str);

[header, data] = CLW_load(fullfile(load_path, fname));


% plot tap force and tapping onset
time = [0 : header.datasize(6)-1] * header.xstep;

ep_mask = strcmp({header.events.code}, num2str(cond_tap_onset));
ep_c10 = find(ep_mask);

ch = 8; % channel to plot, 'Fz'

% shorter segment
time_start = 9.25;
time_duration = n_cycle_plot * par.cycle_dur;
time_idx = dsearchn(time', [time_start time_start+time_duration]');

% find tap-onset time point
tap_signal = squeeze(data(ep_c10(1), ch, :, :, :, time_idx(1):time_idx(2)));

idxl = tap_signal>=0.1;
idxl(1) = 0;
idx = find(idxl);
yest = tap_signal(idx-1)<0.1; 
onset_time_point = idx(yest);  

f = figure('Position', [148 661 1206 205], 'color', 'white');
pnl = panel(f);
pnl.pack('v', [70, 30]);

ax = pnl(1).select();
plot(time(time_idx(1):time_idx(2)), ...
    squeeze(data(ep_c10(1), ch, :, :, :, time_idx(1):time_idx(2))), ...
    'linew', 2);
hold(ax, 'on');
box(ax, 'off');
axis(ax, 'off');
ax.XTick = []; 
ax.YTick = []; 

ax = pnl(2).select();
time_stim = [0 : length(stim.s{cond_tap_onset})-1] * (1/stim.fs);
time_stim_idx = dsearchn(time_stim', [time_start time_start+time_duration]');
plot(time_stim(time_stim_idx(1):time_stim_idx(2)), ...
    stim.s{cond_tap_onset}(time_stim_idx(1):time_stim_idx(2)), ...
    'Color', [128 128 128]/255 , 'linew', 2); %[0.9098 0.4824 0.0549]

box(ax, 'off');
axis(ax, 'off');
ax.XTick = []; 
ax.YTick = []; 

pnl.de.margin = [0, 0, 0, 8];
pnl.margin = [20, 2, 2, 15];


if fig_save
    fname_save = sprintf('TAP-ONSET_extraction_no-onsets_%s_cond-%d', sub_str, cond_tap_onset);
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_folder, fname_save), '.png']);
end


%% plot ITI by IOI
% load tap ratios for each subject from results stored in summary folder

load_path = fullfile(par.summary_path, ...
    sprintf('response-%s_stim-%s', response, feat_name_stim));
fname = sprintf('response-%s_stim-%s_data.mat', response, feat_name_stim);

res_tap_onsets = load(fullfile(load_path, fname));

for i_sub = 1:length(res_tap_onsets.res)
    all_tap_ratios(:, i_sub) = res_tap_onsets.res(i_sub).tap_ratios(:, 1);

end


% plotting
f = figure('pos',[25 485 182 151], 'color', 'white');

pnl = panel(f);
ax = pnl.select();

% color = [49, 97, 127] / 255;
color_group = [136, 34, 85] / 255;
color_ind = [204, 102, 119] / 255;


tested_ratios = stim.ratios(:,1);
mu_tap = mean(all_tap_ratios, 2);
sem_tap = std(all_tap_ratios, [], 2) ./ sqrt(size(all_tap_ratios, 2));

x = [tested_ratios', fliplr(tested_ratios')];
y = [(mu_tap + sem_tap)' , fliplr((mu_tap - sem_tap)')];

lims = [min(tested_ratios) - 0.02, max(tested_ratios) + 0.02];
plot(lims, lims, '--', 'Color', [0.5 0.5 0.5], 'LineWidth', 1.5);
hold on;
fill(x, y, color_ind, 'FaceAlpha', .5, 'EdgeColor', 'none');
hold on;

plot(tested_ratios, mu_tap, 'LineWidth', 2.5, 'Color', color_group);
set(gca, 'YDir','reverse');

ax.XLim = [min(stim.ratios(:, 1)) - 0.02, max(stim.ratios(:, 1)) + 0.02]; 
ax.YLim = [min(stim.ratios(:, 1)) - 0.02, max(stim.ratios(:, 1)) + 0.02]; 

ax.YAxis.LineWidth = 2; 
ax.XAxis.LineWidth = 2; 
ax.XAxis.TickLength = [0,0]; 
ax.YAxis.TickLength = [0,0]; 
ax.XTick = round([min(tested_ratios), max(tested_ratios)], 2); 
ax.YTick = round([min(tested_ratios), max(tested_ratios)], 2); 

pnl.de.margin = 0;
axis(ax, 'square')

if fig_save
    fname_save = sprintf('TAP-ONSET_%s_ratio-2D-plot', response);
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_folder, fname_save), '.png']);
end


%% plot ITI ratio results: RSM and correlation coefficients
% plot overlay best categorical model, all participants
f = plot_rdms_overlay({res_tap_onsets.res.best_model}, 'alpha', 0.1);
title('best m. tap onsets');


if fig_save
    fname_save = 'TAP-ONSET_ratio_overlay-best-model';
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_folder, fname_save), '.png']);
end


%% correlation coefficients ITI ratio
f = figure('pos', [827 508 153 189], 'color', 'white');

pnl = panel(f);
ax = pnl.select();

plot_r_bar([res_tap_onsets.res.r], [res_tap_onsets.res.p_perm], ...
    'ax', ax, ...
    'c_group', [136, 34, 85] / 255, 'c_ind', [204, 102, 119] / 255);

ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
ax.XAxisLocation = 'origin'; 
ax.YTick = [0, 1]; 
ax.FontSize = 14; 
ax.FontName = 'Arial';
ylabel(ax, 'ρ', 'FontName', 'Arial', 'FontSize', 14);

if fig_save
    fname_save = 'TAP-ONSET_ratio_rho_group';
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_folder, fname_save), '.png']);
end