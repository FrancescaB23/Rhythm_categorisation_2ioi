% script to identify consecutive significant harmonics of the continuous 
% tapping force signal at the group level

clc
clear
close all


%% parameters
par = get_par('experiment', '2ioi_conga');
subjects = [num2cell([2:7 9:15 17:20])];

z_thr = norminv(1 - 1e-3); % zsnr significance threshold


%%
data_all = cell(1, length(subjects));

for i_sub = 1:length(subjects)

    sub = subjects{i_sub};
    sub_str = sub_num2str(sub);

    % load data
    fpath = fullfile(par.deriv_path, 'tapping', sub_str);
    fname = sprintf('sub-%03d_task-TAP_response-force.lw6', sub);
    d = dir(fullfile(fpath, fname));

    [header_tap, data_tap] = CLW_load(fullfile(d.folder, d.name));

    % keep only Fz channel
    [header_tap, data_tap] = ...
        RLW_arrange_channels(header_tap, data_tap, {'Fz'});

    % append the data
    data_all{i_sub} = data_tap;
end


%%
% average across all subjects and conditions
data_GA = mean(cat(1, data_all{:}), 1);
header_GA = header_tap;
header_GA.name = sprintf('sub-grand%s', header_tap.name(8:end));
header_GA.datasize = size(data_GA);


% now that I have one signal in the time domain with averaged participants
% and conditions, I compute magnitude FFT
[header_mX, data_mX] = RLW_FFT(header_GA, data_GA);

fmax = (1/header_GA.xstep) / 2;
[header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax);
data_mX(:, :, :, :, :, 1) = 0;
header_mX.name = sprintf('%s_snr-none_mX', header_mX.name);


% compute z-score
[header_GA_z, data_GA_z] = RLW_SNR(header_mX, data_mX, ...
    'xstart', par.snr_bins(1), 'xend', par.snr_bins(2), ...
    'operation', 'zscore');


%% plot zSNR

% highest harmonic
max_harm = floor((fmax-1)/(1/par.cycle_dur));

harm = [1 : max_harm];
frex = harm * (1/par.cycle_dur);
freq = [0 : header_GA_z.datasize(end)-1] * header_GA_z.xstep;

frex_idx = dsearchn(freq', frex');
max_x_plot = dsearchn(freq', 128);

f = figure();
plot(freq(1:max_x_plot), squeeze(data_GA_z(1:max_x_plot)));
hold on;
plot(freq(frex_idx), squeeze(data_GA_z(frex_idx)), 'o');
hold on;
plot([0, 128], [z_thr, z_thr], ...
    '--', 'Color', [0.5, 0.5, 0.5], 'linew', 2);

txt = text(par.filt_cut_high, max(squeeze(data_GA_z(frex_idx))), sprintf('zSNR thr = %.2f', z_thr));
txt.HorizontalAlignment = 'right';
txt.FontSize = 12;

% add markers highlighting harmonics position and color legend for
% significance according to different thresholds for z
h_thr_309 = find(squeeze(data_GA_z(frex_idx))>(norminv(1 - 1e-3)))';
h_thr_232 = find(squeeze(data_GA_z(frex_idx))>(norminv(1 - 1e-2)))';
h_thr_164 = find(squeeze(data_GA_z(frex_idx))>(norminv(1 - .05)))';

pointer_height = floor(max(squeeze(data_GA_z(1:max_x_plot))))+1;
if isempty(h_thr_309)
    h1 = plot(nan, nan, 'v', 'MarkerSize', 10, ...
        'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
else
    h1 = plot(freq(frex_idx(h_thr_309)), ...
        (pointer_height+3)*ones(1,length(h_thr_309)), 'v', ...
        'MarkerSize', 10, 'MarkerEdgeColor', 'r', 'MarkerFaceColor', 'r');
end

if isempty(h_thr_232)
    h2 = plot(nan, nan, 'v', 'MarkerSize', 10, ...
        'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
else
    h2 = plot(freq(frex_idx(h_thr_232)), ...
        (pointer_height+2)*ones(1,length(h_thr_232)), 'v', ...
        'MarkerSize', 10, 'MarkerEdgeColor', 'b', 'MarkerFaceColor', 'b');
end

if isempty(h_thr_164)
    h3 = plot(nan, nan, 'v', 'MarkerSize', 10, ...
        'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
else
    h3 = plot(freq(frex_idx(h_thr_164)), ...
        (pointer_height+1)*ones(1,length(h_thr_164)), 'v', ...
        'MarkerSize', 10, 'MarkerEdgeColor', 'k', 'MarkerFaceColor', 'k');
end

h4 = plot(freq(frex_idx), pointer_height*ones(1,length(frex_idx)), ...
    'v', 'MarkerSize', 10, 'MarkerEdgeColor', 'k');

legend([h1, h2, h3, h4], ...
    {'zSNR>3.09', 'zSNR>2.32', 'zSNR>1.64', 'all harm'}, ...
    'Location', 'best');
legend('Box', 'off');

box off;

title(sprintf('%s GA_%d-subjects %s TAP', ...
    par.experiment, length(subjects)), ...
    'Interpreter', 'none');

zSNR_table = table(harm', frex', squeeze(data_GA_z(frex_idx)));
zSNR_table.Properties.VariableNames = {'harm', 'freq', 'zSNR'};


%%
T = zSNR_table;
T.freq = round(T.freq, 2);
T.zSNR = round(T.zSNR, 2);

find(squeeze(data_GA_z(frex_idx))>z_thr)'
