% This script computes zSNR at the electrodes and frequencies of interest
% This information can be used to select subjects on the base of zSNR
% It also computes a figure for supplementary material showing the chunk
% spectrum across frequency of interest for selected / not selected
% participants

clc
clear all
close all


%% parameters
par = get_par('experiment', 'newborn'); 

subjects = [num2cell([1:3 5:14 16:25 27:29 31:42])];

harmonics = [1:6]; 

roi_name = 'front'; 

z_threshold = norminv(1 - .1); % 1.28 p =.1

fig_save = true;


%%

frex = 1 / par.cycle_dur * harmonics; % frequencies of interest

% initialise max and min values for plotting
% so that these parameters are updated in the loop
max_snr_plot = 0;
min_snr_plot = 0;
max_zbl_plot = 0;
min_zbl_plot = 0;
max_sbl_plot = 0;
min_sbl_plot = 0;


%%

for i_sub = 1:length(subjects)


    % load data
    sub = subjects{i_sub}; 
    sub_str = sub_num2str(sub);

    fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
    fname =  sprintf('sub-%03d_task-EEG_filt_ep_AB-100_chanint_rr-m.lw6', sub); 

    d = dir(fullfile(fpath, fname));
    [header, data] = CLW_load(fullfile(d.folder, d.name));

    
    % select and average channels of interest
    if strcmp(roi_name, 'front')

        roi_chans = {'E19' 'E11' 'E4' 'E20' 'E12' 'E5' 'E118' ...
            'E30' 'E13' 'E6' 'E112' 'E105' 'E7' 'E106' 'E29' ...
            'E111' 'E24' 'E124' 'E28' 'E117'};
    else
        error('roi name %s not recognised', roi_name);
    end

    [header, data] = RLW_arrange_channels(header, data, roi_chans);

    [header, data] = RLW_pool_channels(header, data, roi_chans, ...
        'mixed_channel_label', roi_name, 'keep_original_channels', false);


    % average conditions
    data_avg = mean(data, 1);
    header_avg = header;
    header_avg.datasize = size(data_avg);


    % compute amplitude spectrum
    fmax = (1/header_avg.xstep)/2;
    [header_mX, data_mX] = RLW_FFT(header_avg, data_avg);
    [header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax); 
    data_mX(:, :, :, :, :, 1) = 0; 
    header_mX.name = sprintf('%s_snr-none_mX', sub_str); 


    % calculate zSNR
    freq = [0 : header_mX.datasize(end)-1] * header_mX.xstep;   
    
    zSNR_mean(i_sub) = get_z_snr(data_mX, freq, frex, ...
        par.snr_bins(1), par.snr_bins(2)); 



    %% baseline operations on summed chunks around freq of interest
    % for plotting, SI

    [snr_f, snr_chunk_plot, zbl_f, zbl_chunk_plot, ...
        sbl_f, sbl_chunk_plot] = plot_chunk_baseline_operations(...
        data_mX, freq, frex, par.snr_bins(1), par.snr_bins(2));

    % store output for plotting
    bl_output{i_sub}.sub = sub_str;
    bl_output{i_sub}.snr_f = snr_f;
    bl_output{i_sub}.snr_chunk_plot = snr_chunk_plot;
    bl_output{i_sub}.zbl_f = zbl_f;
    bl_output{i_sub}.zbl_chunk_plot = zbl_chunk_plot;
    bl_output{i_sub}.sbl_f = sbl_f;
    bl_output{i_sub}.sbl_chunk_plot = sbl_chunk_plot;

    % update min and max values for plotting
    max_snr_plot = max(max_snr_plot, max(snr_chunk_plot));
    min_snr_plot = min(min_snr_plot, min(snr_chunk_plot));
    max_zbl_plot = max(max_zbl_plot, max(zbl_chunk_plot));
    min_zbl_plot = min(min_zbl_plot, min(zbl_chunk_plot));
    max_sbl_plot = max(max_sbl_plot, max(sbl_chunk_plot));
    min_sbl_plot = min(min_sbl_plot, min(sbl_chunk_plot));

end
    

%% display participant code for participant with zSNR > threshold

zSNR_mean

subjects(zSNR_mean > z_threshold) 

sub_for_selection = cell2mat(subjects(zSNR_mean > z_threshold));
disp(sprintf('subjects with zSNR > %.2f', z_threshold));
fprintf('%d ', sub_for_selection);


T = table(cell2mat(subjects)', zSNR_mean', 'VariableNames', {'sub', 'zSNR'})


%% supplementary figure: zSNR at the individual level
% overimposed across selected (left) and not selected (right) participants

f = figure('Position', [252 580 554 290]); % [253 577 783 293]

xtick_pos = [1, par.snr_bins(2)+1, par.snr_bins(2)*2+1];
xtick_labels = {round(-header_mX.xstep*par.snr_bins(2), 2), 'f', ...
    round(header_mX.xstep*par.snr_bins(2), 2)};

ytick_pos = [-2, z_threshold, 7];
ytick_labels = {ytick_pos(1), round(z_threshold, 2), ytick_pos(end)};

% plot zSNR of participants exceeding threshold at the frequency of
% interest
subplot(1, 2, 1)
hold on;
for i_sub = 1:size(bl_output, 2)

    if bl_output{i_sub}.zbl_f > z_threshold
        plot(bl_output{i_sub}.zbl_chunk_plot, 'LineWidth', 2, 'color', 'b');
    end
end

hold off;
yline(z_threshold, '--k', 'LineWidth', 2);
ylim([min_zbl_plot-0.2 max_zbl_plot+0.2]);
xlim([1 par.snr_bins(2)*2+1]);
ax = gca;
set(ax, 'FontName', 'Arial', 'Fontsize', 12);
xticks(xtick_pos);
xticklabels(xtick_labels);
yticks(ytick_pos);
yticklabels(ytick_labels);
ylabel('z-score');
xlabel('Frequency [Hz]');

% plot zSNR of participants not exceeding threshold at the frequency of
% interest
subplot(1, 2, 2)
hold on;
for i_sub = 1:size(bl_output, 2)

    if ~(bl_output{i_sub}.zbl_f > z_threshold)
        plot(bl_output{i_sub}.zbl_chunk_plot, 'LineWidth', 2, 'color', [.7 .7 .7]);
    end

end

hold off;
yline(z_threshold, '--k', 'LineWidth', 2);
ylim([min_zbl_plot-0.2 max_zbl_plot+0.2]);
xlim([1 par.snr_bins(2)*2+1]);
ax = gca;
set(ax, 'FontName', 'Arial', 'Fontsize', 12);
xticks(xtick_pos);
xticklabels(xtick_labels);
yticks(ytick_pos);
yticklabels(ytick_labels);
ylabel('z-score');
xlabel('Frequency [Hz]');

if fig_save
    folder_save = fullfile(par.deriv_path, 'figures', 'supplementary_material');
    mkdir(folder_save);
    fname_save = sprintf('zSNR_newborn_individual_th-%.2f', z_threshold);
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(folder_save, fname_save), '.png']);
end

