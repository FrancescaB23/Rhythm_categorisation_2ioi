% plot topoplots newborns and adults: specify group names: as the procedure
% is the same regardless of group, like this the same operations are done
% on both datasets without repeating same command lines

clc
clear
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

fig_save = true;


%%
figure_path = fullfile(par_adults.deriv_path, 'figures', 'topoplots');
if ~exist(figure_path)
    mkdir(figure_path);
end

harm_str = vec_to_str(harm_selections, 'format','%d', 'sep','-'); 


%%
data_to_topoplot = {'newborns', 'adults'}; % 'adults' or 'newborns'

for i = 1:2
    % load time domain data, grand average
    switch data_to_topoplot{i}

        case 'adults'

            load_path = fullfile(par_adults.deriv_path, ...
                'preproc', 'sub-grand');
            fname = 'sub-grand_task-EEG_filt_ep_chanint_icfilt_rr-m';
            [header, data] = CLW_load(fullfile(load_path, fname));

        case 'newborns'

            load_path = fullfile(par_newborns.deriv_path, ...
                'preproc', sprintf('sub-grand-%s', ga_newborns));
            fname = sprintf('sub-grand-%s_task-EEG_filt_ep_AB-100_chanint_rr-m', ga_newborns);
            [header, data] = CLW_load(fullfile(load_path, fname));

    end

    % average across conditions
    data = mean(data, 1);
    header.datasize = size(data);

    % compute magnitude spectrum
    fmax = (1/header.xstep) / 2;
    [header_mX, data_mX] = RLW_FFT(header, data);
    [header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax);
    data_mX(:, :, :, :, :, 1) = 0;

    % perform baseline correction
    [header_mX_snr, data_mX_snr] = RLW_SNR(header_mX, data_mX, ...
        'xstart', par_adults.snr_bins(1), 'xend', par_adults.snr_bins(2));

    % sum responses at harmonics of interest
    frex = 1/par_adults.cycle_dur * harm_selections;
    freq = [0 : header_mX_snr.datasize(end)-1] * header_mX_snr.xstep;
    frex_idx = dsearchn(freq', frex');
    data_topoplot = sum(data_mX_snr(:,:,:,:,:,frex_idx), 6);

    % topoplot
    fig = figure('Position', [332 760 143 106])
    fig = CLW_topoplot_vector(header_mX_snr, data_topoplot, ...
        'maplimits', [0 max(data_topoplot)], ...
        'shading', 'interp', 'whitebk', 'on');
    colorbar;
    title(data_to_topoplot{i});

    if fig_save
        print(fig, '-dpng', '-painters', '-r600', fullfile(figure_path, ...
            sprintf('topoplot_%s-%s_h-%s_allcond.png', ...
            data_to_topoplot{i}, feat_resp, harm_str)));
    end

end
