% script to determine the number of significant harmonics (eeg signal)

clc
clear
close all


%% parameters
par = get_par('experiment', '2ioi_conga'); % '2ioi_conga' 'newborn'

switch par.experiment
    case '2ioi_conga'
        subjects = [num2cell([2:7 9:15 17:20])];
    case 'newborn'
        subjects = [num2cell([1:3 5:14 16:25 27:29 31:42])]; % all 38 selected subjects
end


% highest harmonics - dependend on low pass filter parameters
max_harm = floor((par.filt_cut_high-0.5)/(1/par.cycle_dur));

z_thr = norminv(1 - 1e-3); % zsnr significance threshold


%% load individual data

data_all = cell(1, length(subjects));

for i_sub = 1:length(subjects)

    sub = subjects{i_sub};
    sub_str = sub_num2str(sub);


    switch par.experiment

        case 'newborn'
            fpath = fullfile(par.deriv_path, 'preproc', sub_str);
            fname =  sprintf('sub-%03d_task-EEG_filt_ep_AB-100_chanint_rr-m.lw6', sub);

        case '2ioi_conga'
            fpath = fullfile(par.deriv_path, 'preproc', sub_str);
            fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*_rr-m.lw6', sub);
    end

    d = dir(fullfile(fpath, fname));

    [header, data] = CLW_load(fullfile(d.folder, d.name));


    % keep subset of channels for determination of significant harmonics 
    % (depending on group) 
    if strcmp(par.experiment, '2ioi_conga') % 64 channels (all excluding mastoids)

        roi_chans = {header.chanlocs(1:64).labels};
        roi_name = '64chans';

    elseif strcmp(par.experiment, 'newborn') % front channels

        roi_chans = {'E19' 'E11' 'E4' 'E20' 'E12' 'E5' 'E118' ...
            'E30' 'E13' 'E6' 'E112' 'E105' 'E7' 'E106' 'E29' ...
            'E111' 'E24' 'E124' 'E28' 'E117'};
        roi_name = 'front';
    end

    [header, data] = RLW_arrange_channels(header, data, roi_chans);

    % append the data
    data_all{i_sub} = data;

end


%%
% average across all subjects and conditions
data_GA = mean(cat(1, data_all{:}), 1);
header_GA = header;
header_GA.name = sprintf('sub-grand%s', header.name(8:end));
header_GA.datasize = size(data_GA);


% average across all channels
[header_GA, data_GA] = RLW_pool_channels(header_GA, data_GA, roi_chans, ...
    'mixed_channel_label', roi_name, ...
    'keep_original_channels', false);


% now that I have one signal in the time domain with averaged participants,
% conditions and electrodes, I compute magnitude FFT
[header_mX, data_mX] = RLW_FFT(header_GA, data_GA); 

fmax = (1/header_GA.xstep) / 2;
[header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax); 

data_mX(:, :, :, :, :, 1) = 0; 

header_mX.name = sprintf('%s_snr-none_mX', header_mX.name); 


% compute z-score
[header_GA_z, data_GA_z] = RLW_SNR(header_mX, data_mX, ...
    'xstart', par.snr_bins(1), 'xend', par.snr_bins(2), ...
    'operation', 'zscore');


%% identify significant harmonics
harm = [1 : max_harm];
frex = harm * (1/par.cycle_dur);
freq = [0 : header_GA_z.datasize(end)-1] * header_GA_z.xstep;

frex_idx = dsearchn(freq', frex');
max_x_plot = dsearchn(freq', par.filt_cut_high);

f = figure();
plot(freq(1:max_x_plot), squeeze(data_GA_z(1:max_x_plot)));
hold on;
plot(freq(frex_idx), squeeze(data_GA_z(frex_idx)), 'o');
hold on;
plot([0, par.filt_cut_high], [z_thr, z_thr], ...
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

title(sprintf('%s GA_%d-subjects rr-m chans-%s', ...
    par.experiment, length(subjects), roi_name), ...
    'Interpreter', 'none');

zSNR_table = table(harm', frex', squeeze(data_GA_z(frex_idx)));
zSNR_table.Properties.VariableNames = {'harm', 'freq', 'zSNR'};