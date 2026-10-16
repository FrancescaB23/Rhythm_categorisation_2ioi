function rsa_lf_time(sub, par_rsa, varargin)


% function for rsa analysis of frequency domain responses loading time
% domain data

parser = inputParser(); 

addParameter(parser, 'feat_resp', 'X'); 
addParameter(parser, 'feat_stim', 'X'); 
addParameter(parser, 'rerun_stats', true); 
addParameter(parser, 'rerun_figs', true); 
addParameter(parser, 'par', []); 
addParameter(parser, 'flip_rdm', true); 
addParameter(parser, 'reref', []); 


parse(parser, varargin{:}); 

feat_name_resp = parser.Results.feat_resp; 
feat_name_stim = parser.Results.feat_stim; 
rerun_stats = parser.Results.rerun_stats; 
rerun_figs = parser.Results.rerun_figs; 
par = parser.Results.par; 
flip_rdm = parser.Results.flip_rdm; 

% if parameter structure hasn't been explicitly passed, load a new one from
% the file - this allows the possibility for the user to pass an adjusted
% par structure intead of the default one
if isempty(par)
    par = get_par(); 
end


%% 
% frequencies of interest (par_rsa.harmonics of rhythm cycle repetition rate) 
harm_str = vec_to_str(par_rsa.harmonics, 'format','%d', 'sep','-'); 

frex = 1 / par.cycle_dur * par_rsa.harmonics; 


%% input directories 

sub_str = sub_num2str(sub); 
load_path = fullfile(par.deriv_path,'preproc', sub_str); 


%% output directories

if par_rsa.feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end

subdirs = fullfile(sprintf('response-lf_stim-%s', feat_name_stim), ...
                     sprintf('harm-%s', harm_str), ...
                     sprintf('snr-%s', snr_str), ...
                     sprintf('roi-%s_avgChans-%s', par_rsa.roi_name, jsonencode(par_rsa.avg_chans)) ...
                     ); 

save_path = fullfile(par.results_path, subdirs, sub_str); 

mkdir(save_path); 

    
%% open csv file 
fname_tbl = sprintf('response-lf_feat-%s_harm-%s_snr-%s_roi-%s_avgChans-%s_cat-2-3-4-5.csv', ...
    feat_name_resp, ...
    harm_str, ...
    snr_str, ...
    par_rsa.roi_name, ...
    jsonencode(par_rsa.avg_chans) ...
    );

if exist(fullfile(par.results_path, subdirs, fname_tbl))
    tbl_res = readtable(fullfile(par.results_path, subdirs, fname_tbl));
else
    col_names = {'sub', 'n_categ_fit', 'bound1', 'bound2', 'bound3', ...
        'bound4', 'r', 'p', 'p_adj', ...
        'z_snr', 'z_snr_sd', 'n_models_fit', 'n_perm', 'chans'};
    tbl_res = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);
end


%% load stimulus info 

load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'));

if strcmp(feat_name_stim, 'cochlear') % cochlear model as regressor

    % load cochlear fft responses
    fpath_cochlear = fullfile(['/Users/barbero/Documents/Projects/' ...
        'Rhototype_2ioi_newborns/adults/data/derivatives'], ...
        'fft', 'response-urear-an');
    fname_cochlear = sprintf('response-urear-an_snr-none_%s', feat_name_resp);
    [header_cochlear, data_cochlear] = CLW_load(fullfile(fpath_cochlear, fname_cochlear));

    % get frequency vector - specific for cochlear responses
    freq_cochlear = [0: header_cochlear.datasize(end)-1] * header_cochlear.xstep;

    % get cochlear rdm
    feat_for_rdm_cochlear = get_feat_from_fft(header_cochlear, data_cochlear, freq_cochlear, frex);

    rdm_stim = get_rdm_from_feat(feat_for_rdm_cochlear, 'method', 'pearson');

    figure()
    imagesc(rdm_stim)
    title(sprintf('cochlear %s', feat_name_resp))

else

    rdm_stim = rsa_stim_lf(par, frex, 'feat', feat_name_stim);

    figure()
    imagesc(rdm_stim);
    title(sprintf('%s', feat_name_stim));
end


%% load time domain data

if strcmp(par.experiment, 'newborn')
    fname = sprintf('%s_task-EEG_filt_ep_AB-100_chanint_rr-m', sub_str);

    [header, data] = CLW_load(fullfile(load_path, fname));

else 

    fname = sprintf('%s_task-EEG_filt_ep_chanint_icfilt*rr-m.lw6', sub_str);
    d = dir(fullfile(load_path, fname)); 

    [header, data] = CLW_load(fullfile(d.folder, d.name));
end


%% compute magnitude FFT with snr correction for topoplot

fmax = (1/header.xstep) / 2;

[header_mX_allch, data_mX_allch] = RLW_FFT(header, data);

[header_mX_allch, data_mX_allch] = crop_lw_data(header_mX_allch, data_mX_allch, fmax);

data_mX_allch(:, :, :, :, :, 1) = 0;

header_mX_allch.name = sprintf('%s_snr-none_mX', sub_str);


% subtract surrounding bins
[header_mX_allch_snr, data_mX_allch_snr] = RLW_SNR( ...
    header_mX_allch, data_mX_allch, ...
    'xstart', par.snr_bins(1), 'xend', par.snr_bins(2));

header_mX_allch_snr.name = sprintf('%s_snr-%d-%d_mX', sub_str, ...
    par.snr_bins(1), par.snr_bins(2));


% save a copy of all chanlocs for topoplot
chanlocs_all_chans = header_mX_allch.chanlocs;


% get amplitude channels for topoplots and zsnr all chans. zsnr will be
% recomputed later if the roi considers average channels
freq = [0 : header_mX_allch.datasize(end)-1] * header_mX_allch.xstep;   

z_snr_all_chan = get_z_snr(data_mX_allch, freq, frex, ...
    par.snr_bins(1), par.snr_bins(2)); 

amp_sum_all_chan = get_amp_summary(data_mX_allch_snr, ...
    freq, frex, 'method', 'sum'); 


%% channel selection and average
if strcmp(par_rsa.roi_name, 'front')

    if strcmp(par.experiment, 'newborn')

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


    else

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
end

roi_chans_idx = get_chan_idx(header, roi_chans); 

roi_name_for_tbl = {'front'}; 

roi_chans


[header, data] = RLW_arrange_channels(header, data, roi_chans);

if par_rsa.avg_chans

    [header, data] = RLW_pool_channels(header, data, roi_chans, ...
        'mixed_channel_label', roi_name_for_tbl{:}, ...
        'keep_original_channels', false);
end


%% compute amplitude and complex FFT

[header_mX, data_mX] = RLW_FFT(header, data);

[header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax);

data_mX(:, :, :, :, :, 1) = 0;

header_mX.name = sprintf('%s_snr-none_mX', sub_str);


[header_X, data_X] = RLW_FFT(header, data, 'output', 'complex');

[header_X, data_X] = crop_lw_data(header_X, data_X, fmax);

data_X(:, :, :, :, :, 1) = 0;

header_X.name = sprintf('%s_snr-none_X', sub_str);


% noise subtracted spectra
if par_rsa.feat_snr

    % magnitude spectra
    [header_mX_snr, data_mX_snr] = RLW_SNR(header_mX, data_mX,...
        'xstart', par.snr_bins(1), 'xend', par.snr_bins(2));

    header_mX_snr.name = sprintf('%s_snr-%d-%d_mX', ...
        sub_str, par.snr_bins(1), par.snr_bins(2));


    % complex spectra
    header_X_snr = header_X;

    [data_X_snr] = subtract_noise_bins_complex(...
        data_X, par.snr_bins(1), par.snr_bins(2));

    header_X_snr.name = sprintf('%s_snr-%d-%d_X', ...
        sub_str, par.snr_bins(1), par.snr_bins(2));

end


%% re-estimate zSNR
z_snr_roi = get_z_snr(data_mX, freq, frex, par.snr_bins(1), par.snr_bins(2)); 
z_snr_roi = mean(z_snr_roi, 2); 


%% get RDM

switch feat_name_resp

    case 'mX'

        if par_rsa.feat_snr
            feat_for_rdm = get_feat_from_fft(header_mX_snr, data_mX_snr, freq, frex);
        else
            feat_for_rdm = get_feat_from_fft(header_mX, data_mX, freq, frex);
        end

    case 'X'

        if par_rsa.feat_snr
            feat_for_rdm = get_feat_from_fft(header_X_snr, data_X_snr, freq, frex);
        else
            feat_for_rdm = get_feat_from_fft(header_X, data_X, freq, frex);
        end
end

rdm = get_rdm_from_feat(feat_for_rdm, 'method', 'pearson'); 


%% plotting 

% plot topo 
fname = sprintf('%s_feat-ampSum_topo', sub_str); 

if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs


    cmap = parula(); 

    % plot overall response magnitude disstribution across scalp 
    f = plot_topoplots(amp_sum_all_chan, chanlocs_all_chans, ...
                  'lab', ...
                  'amp-sum', ...
                  'mark_chan_idx', roi_chans_idx, ...
                  'cmap', cmap); 

    print(f, '-dpng', '-painters', '-r90', fullfile(save_path, fname));

    close all
end

% plot magnitude spectra 
fname = sprintf('%s_mX', sub_str); 

if par_rsa.feat_snr
    header_mX_to_plot = header_mX_snr; 
    data_mX_to_plot = data_mX_snr; 
else
    header_mX_to_plot = header_mX; 
    data_mX_to_plot = data_mX; 
end

if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs

    f = plot_mX_lw(header_mX_to_plot, data_mX_to_plot, ...
                   [1 : header_mX_to_plot.datasize(2)], ...
                   frex, ...
                   'z_snr', z_snr_roi, ...
                   'fmax', max([frex, 15])); 

    print(f, '-dpng', '-painters', '-r200', fullfile(save_path, fname))

end

data_to_save = []; 
data_to_save.header_mX_to_plot = header_mX_to_plot; 
data_to_save.data_mX_to_plot = data_mX_to_plot; 
data_to_save.frex = frex; 
data_to_save.z_snr = z_snr_roi; 

save(fullfile(save_path, [fname, '.mat']), '-struct', 'data_to_save'); 


% plot RDM 
fname = sprintf('%s_feat-%s_method-pearson_rdm', sub_str, feat_name_resp); 

if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs

    f = plot_rdm(rdm, 'ratios', stim.ratios, 'flip', flip_rdm); 

    print(f, '-dpng', '-painters', '-r200', fullfile(save_path, fname))

    close all

end

save(fullfile(save_path, [fname, '.mat']), 'rdm'); 



%% 2, 3, 4, or 5 categories models

n_categ = [2, 3, 4, 5]; 

fname = sprintf('%s_feat-%s_method-pearson_nCateg-2-3-4-5', sub_str, feat_name_resp); 

fname_mat = fullfile(save_path, [fname, '_perm.mat']); 

if isfile(fname_mat) & ~rerun_stats
    
    warning('loading permutation data from file `%s`', fname); 
    
    perm_output = load(fname_mat); 
    
    r_obs = perm_output.r_obs; 
    p_perm = perm_output.p_perm; 
    r_null_dist = perm_output.r_null_dist; 
    best_model = perm_output.best_model; 
    all_models_info = perm_output.all_models_info; 
    
else
    
    [r_obs, p_perm, r_null_dist, best_model, best_model_info, all_models_info] = ...
                          find_test_best_categ_model(...
                                        rdm, n_categ, par.n_perm_rsa_individual, ...
                                        'covariate', rdm_stim, ...
                                        'n_mini_categ_skip', par.n_mini_categ_skip);
                                    
    save(fullfile(save_path, [fname, '_perm.mat']), ...
                  'par', 'par_rsa', 'r_obs', 'p_perm', 'r_null_dist', ...
                  'best_model', 'all_models_info');  
                         
    writetable(all_models_info, fullfile(save_path, ...
               [fname, '_models.csv'])); 

end

% write table 
if (isnumeric(sub) || contains(sub_str, 'ga1out')) && ~contains(sub_str, 'grand')

    if isnumeric(sub)
        row_mask = (tbl_res.sub == sub) & (tbl_res.n_categ_fit == 5); 
    elseif contains(sub, 'ga1out')
        row_mask = strcmp(tbl_res.sub, sub) & (tbl_res.n_categ_fit == 5); 
    end 
    tbl_res(row_mask, :) = []; 

    new_row = [...
        {sub}, ...
        5, ...
        all_models_info{1, 'bound1'}, ...
        all_models_info{1, 'bound2'}, ...
        all_models_info{1, 'bound3'}, ...
        all_models_info{1, 'bound4'}, ...
        all_models_info{1, 'r'}, ...
        all_models_info{1, 'p'}, ...
        all_models_info{1, 'p_adj'}, ...
        mean(z_snr_roi), ...
        std(z_snr_roi), ...
        size(all_models_info, 1), ...
        par.n_perm_rsa_individual, ...
        roi_name_for_tbl, ...
        ]; 

    tbl_res = [tbl_res; new_row]; 
end

% plot 
if ~isfile(fullfile(save_path, [fname, '_perm.png'])) || rerun_figs

    f = plot_best_model(rdm, best_model, r_obs, p_perm, r_null_dist, ...
                        'z_snr', z_snr_roi, 'ratios', stim.ratios, 'flip', flip_rdm);

    print(f, '-dpng', '-painters', '-r100', ...
          fullfile(save_path, [fname, '_perm']))

    close all
end


%% 
if isnumeric(sub)
    writetable(tbl_res, fullfile(par.results_path, subdirs, fname_tbl));
end
end
