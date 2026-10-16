function rsa_tap_force(sub, par_rsa, varargin)

% function for fRSA analysis of frequency domain tapping force responses
% loading time domain data

parser = inputParser();

addParameter(parser, 'feat_resp', 'X');
addParameter(parser, 'feat_stim', 'X');
addParameter(parser, 'rerun_stats', true);
addParameter(parser, 'rerun_figs', true);
addParameter(parser, 'par', []);
addParameter(parser, 'flip_rdm', true);

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
load_path = fullfile(par.deriv_path, 'tapping', sub_str);


%% output directories

if par_rsa.feat_snr
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2));
else
    snr_str = 'none';
end

data_path = fullfile(par.results_path, ...
    sprintf('response-tap-force_stim-%s', feat_name_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str));

save_path = fullfile(data_path, sub_str);

mkdir(save_path);

%% open csv file
fname_tbl = sprintf('response-tap-force_feat-%s_harm-%s_snr-%s_cat-2-3-4-5.csv', ...
    feat_name_resp, ...
    harm_str, ...
    snr_str);

if exist(fullfile(data_path, fname_tbl))
    tbl_res = readtable(fullfile(data_path, fname_tbl));
else
    col_names = {'sub', 'n_categ_fit', 'bound1', 'bound2', 'bound3', ...
        'bound4', 'r', 'p', 'p_adj', 'z_snr', 'z_snr_sd', ...
        'n_models_fit', 'n_perm'};
    tbl_res = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);
end


%% load stimulus info

load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'));
rdm_stim = rsa_stim_lf(par, frex, 'feat', feat_name_stim);


%% load time domain data

fname = sprintf('%s_task-TAP_response-force_avg.lw6', sub_str);
[header_tap, data_tap] = CLW_load(fullfile(load_path, fname));


%% compute FFT
fmax = (1/header_tap.xstep) / 2;

[header_mX, data_mX] = RLW_FFT(header_tap, data_tap);

[header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax);

data_mX(:, :, :, :, :, 1) = 0;

header_mX.name = sprintf('%s_snr-none_mX', sub_str);


[header_X, data_X] = RLW_FFT(header_tap, data_tap, 'output', 'complex');

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


%% get z_snr

freq = [0 : header_mX.datasize(end)-1] * header_mX.xstep;

% get zSNR from raw magnitude spectra
z_snr = get_z_snr(data_mX, freq, frex, par.snr_bins(1), par.snr_bins(2));


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

% plot magmitude spectra
fname = sprintf('%s_mX', sub_str);

if par_rsa.feat_snr
    header_mX_to_plot = header_mX_snr;
    data_mX_to_plot = data_mX_snr;
else
    header_mX_to_plot = header_mX;
    data_mX_to_plot = data_mX;
end

if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs

    f = plot_mX_lw(header_mX_to_plot, data_mX_to_plot, 1, frex, ...
        'z_snr', z_snr, 'fmax', max([frex, 15]));

    print(f, '-dpng', '-painters', '-r600', fullfile(save_path, fname))

end

data_to_save = [];
data_to_save.header_mX_to_plot = header_mX_to_plot;
data_to_save.data_mX_to_plot = data_mX_to_plot;
data_to_save.frex = frex;
data_to_save.z_snr = z_snr;

save(fullfile(save_path, [fname, '.mat']), '-struct', 'data_to_save');


% plot RDM
fname = sprintf('%s_feat-%s_method-pearson_rdm', sub_str, feat_name_resp);

if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs

    f = plot_rdm(rdm, 'ratios', stim.ratios, 'flip', flip_rdm);

    print(f, '-dpng', '-painters', '-r300', fullfile(save_path, fname))
    save(fullfile(save_path, [fname, '.mat']), 'rdm');

    close all
end


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
if isnumeric(sub)
    row_mask = (tbl_res.sub == sub) & (tbl_res.n_categ_fit == 5);
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
        mean(z_snr), ...
        std(z_snr), ...
        size(all_models_info, 1), ...
        par.n_perm_rsa_individual, ...
        ];

    tbl_res = [tbl_res; new_row];
end

if ~isfile(fullfile(save_path, [fname, '_perm.png'])) || rerun_figs

    f = plot_best_model(rdm, best_model, r_obs, p_perm, r_null_dist, ...
        'z_snr', z_snr, 'ratios', stim.ratios, 'flip', flip_rdm);

    print(f, '-dpng', '-painters', '-r300', ...
        fullfile(save_path, [fname, '_perm']))

    close all

end


%%
if isnumeric(sub)
    writetable(tbl_res, fullfile(data_path, fname_tbl));
end
