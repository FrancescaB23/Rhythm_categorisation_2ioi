function rsa_tap_ratios(sub, varargin)

parser = inputParser();

addParameter(parser, 'feat_stim', 'ratio');
addParameter(parser, 'rerun_stats', true);
addParameter(parser, 'rerun_figs', true);
addParameter(parser, 'par', []);
addParameter(parser, 'flip_rdm', true);

parse(parser, varargin{:});

feat_stim = parser.Results.feat_stim;
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


%% input directories
sub_str = sub_num2str(sub);
load_path = fullfile(par.deriv_path, 'tap_onsets', sub_str);


%% output directories
data_path = fullfile(par.results_path, sprintf('response-tap-onsets_stim-%s', feat_stim));
save_path = fullfile(data_path, sub_str);
mkdir(save_path);


%% open csv file
fname_tbl = sprintf('response-tap-onsets_cat-2-3-4-5.csv');

if exist(fullfile(data_path, fname_tbl))
    tbl_res = readtable(fullfile(data_path, fname_tbl));
else
    col_names = {'sub', 'n_categ_fit', 'bound1', 'bound2','bound3', ...
        'bound4', 'r', 'p', 'p_adj', 'n_models_fit', 'n_perm'};
    tbl_res = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);
end


%% load data
fname = sprintf('%s_response-tap-ratios.mat', sub_str);
data = load(fullfile(load_path, fname));


%% load stimulus info

load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'));

if strcmp(feat_stim, 'ratio')

    rdm_stim = rsa_stim_lf(par, [], 'feat', feat_stim);

else

    error('feat_stim = %s not implementeed', feat_stim);
end

%% get RDM
rdm = get_rdm_from_feat(data.tap_ratios', 'method', 'L2');


%%
fname = sprintf('%s_feat-ratios_method-pearson_rdm', sub_str);

if ~isfile(fullfile(save_path, [fname, '.png'])) || rerun_figs

    f = plot_rdm(rdm, 'ratios', stim.ratios, 'flip', flip_rdm);

    print(f, '-dpng', '-painters', '-r300', fullfile(save_path, fname))
    save(fullfile(save_path, [fname, '.mat']), 'rdm');

    close all
end


%% 2, 3, 4, or 5 categories models

n_categ = [2, 3, 4, 5];

fname = sprintf('%s_feat-ratios_method-pearson_nCateg-2-3-4-5', sub_str);

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
        'par', 'r_obs', 'p_perm', 'r_null_dist', ...
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
        size(all_models_info, 1), ...
        par.n_perm_rsa_individual, ...
        ];

    tbl_res = [tbl_res; new_row];
end

% plot
if ~isfile(fullfile(save_path, [fname, '_perm.png'])) || rerun_figs

    f = plot_best_model(rdm, best_model, r_obs, p_perm, r_null_dist, ...
        'ratios', stim.ratios, 'flip', flip_rdm);

    print(f, '-dpng', '-painters', '-r300', ...
        fullfile(save_path, [fname, '_perm']))

    close all

end

%%
if isnumeric(sub)
    writetable(tbl_res, fullfile(data_path, fname_tbl));
end

