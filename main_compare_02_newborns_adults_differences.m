%% Explore differences in representational geometries for EEG in newborns and adults
% Differences are explored in two ways:
% 1) a quantitative analyses testing one measure on the best model from
% another measure
% 2) qualitatively using multidimensional scaling and hierarchical
% clustering

clear
close all
clc


%% parameters
par_adults = get_par('experiment', '2ioi_conga');
par_newborns = get_par('experiment', 'newborn');

ga_newborns = '19'; 

harm_selections = [1:6];
feat_snr = false;

roi_name = 'front';
avg_chans = true;

feat_resp = 'mX'; 

feat_stim = 'ratio';

n_perm = 10000;

fig_save = false;
run_stats = false;


%%
harm_str = vec_to_str(harm_selections, 'format','%d', 'sep','-'); 

if feat_snr
    snr_str = sprintf('%d-%d', par_adults.snr_bins(1), par_adults.snr_bins(2));
else
    snr_str = 'none';
end

save_path = fullfile(par_adults.deriv_path, 'comparisons', ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)));
mkdir(save_path);


%% load results
% load results adults
res_adults = load(fullfile( ...
    par_adults.summary_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf(['response-lf_feat-%s_harm-%s_snr-%s_' ...
    'roi-%s_avgChans-%s_cat-2-3-4-5_data.mat'], ...
    feat_resp, harm_str, snr_str, roi_name, jsonencode(avg_chans))));


% load results newborns
res_newborns = load(fullfile( ...
    par_newborns.results_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf('sub-grand-%s', ga_newborns), ...
    sprintf('sub-grand-%s_feat-%s_method-pearson_nCateg-2-3-4-5_perm.mat', ...
    ga_newborns, feat_resp)));

res_newborns.rdm = load(fullfile( ...
    par_newborns.results_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf('sub-grand-%s', ga_newborns), ...
    sprintf('sub-grand-%s_feat-%s_method-pearson_rdm.mat', ...
    ga_newborns, feat_resp)));


%% save fname
fname_save_base = fullfile(save_path, ...
    sprintf('response-lf_feat-%s_featStim-%s_harm-%s_snr-%s_roi-%s_avgChans-%s_adultsVsNewborns', ...
    feat_resp, ...
    feat_stim, ...
    harm_str, ...
    snr_str, ...
    roi_name, ...
    jsonencode(avg_chans)));


%% Define "winning model" for adults
% Approach: "most popular" theoretical model

x = [res_adults.res.best_model_id];

% Get unique values and the number of their occurrence
[unique_vals, ~, ic] = unique(x);
counts = accumarray(ic, 1);

% sort by counts in descending order
[sorted_counts, sort_idx] = sort(counts, 'descend');
sorted_values = unique_vals(sort_idx)';

disp('model_id  occurrences');
disp([sorted_values, sorted_counts]);

model_id_adults_A = sorted_values(1);


%% Define "winning model" for newborns

r_newborns_all = sortrows(res_newborns.all_models_info, ...
    {'n_categ', 'bound1', 'bound2', 'bound3', 'bound4'}).r;

[max_r_newborns, model_id_newborns] = max(r_newborns_all);


%% stimuli info

load(fullfile(par_adults.deriv_path, 'stimuli', 'stim.mat')); 
rdm_stim = rsa_stim_lf(par_adults, [], 'feat', feat_stim);

figure();
imagesc(rdm_stim); title(sprintf('%s', feat_stim));


%% categorical models
% generate all possible categorical models
n_categ = [2, 3, 4, 5];
n_mini_categ_skip = 1;
[models, all_models_info] = generate_categ_models_5(...
    par_adults.n_cond, n_categ, ...
    'n_mini_categ_skip', n_mini_categ_skip);


%% sanity check: regenerated models order match the order of the results
cols_check = {'n_categ', 'bound1', 'bound2', 'bound3', 'bound4'};

if isequaln(res_adults.res(1).all_models_info{:, cols_check}, all_models_info{:, cols_check})
    disp('model definitions match across adults and re-generated models, continuing the script');
else
    error('model definitions do not match the regenerated ones, check!');
end


%% Test newborns data on adult eeg winning model

if exist([fname_save_base, '_res_n2a.mat']) && ~run_stats
    clear res_n2a
    load([fname_save_base, '_res_n2a.mat']);
else
    res_n2a = [];
    
    [res_n2a.r_obs, res_n2a.p_perm, res_n2a.r_null_dist] = ...
        perm_rdm_corr(...
        res_newborns.rdm.rdm, ...
        models{model_id_adults_A}, ...
        n_perm, ...
        'parfor', true, ...
        'covariate', rdm_stim);
    
    save([fname_save_base, '_res_n2a.mat'], 'res_n2a');
end


f = plot_best_model(res_newborns.rdm.rdm, models{model_id_adults_A}, ...
    res_n2a.r_obs, res_n2a.p_perm, res_n2a.r_null_dist, ...
    'ratios', stim.ratios);
sgtitle('newborns on adult model')
if fig_save
    saveas(f, [fname_save_base, '_res_n2a.png']);
end



%% Test adult data on newborns eeg winning model

% group-level permutation test
if exist([fname_save_base, '_res_a2n_group.mat']) && ~run_stats
    clear res_a2n_group
    load([fname_save_base, '_res_a2n_group.mat']);
else
    res_a2n_group = [];
    [res_a2n_group.r_obs, res_a2n_group.r_null_dist, res_a2n_group.p_perm] = ...
        perm_rdm_corr_group_empirical_model(...
        {res_adults.res.rdm}, ...
        repmat({models{model_id_newborns}}, 1, size({res_adults.res.rdm}, 2)), ...
        n_perm, ...
        'parfor', true, ...
        'covariate', rdm_stim);
    save([fname_save_base, '_res_a2n_group.mat'], 'res_a2n_group');
end

% individual level permutation test with fixed model
if exist([fname_save_base, '_res_a2n_ind.mat']) && ~run_stats
    clear res_a2n_ind
    load([fname_save_base, '_res_a2n_ind.mat']);
else
    res_a2n_ind = [];
    for i_sub=1:length(res_adults.res)
        i_sub
        [res_a2n_ind(i_sub).r_obs, res_a2n_ind(i_sub).p_perm, res_a2n_ind(i_sub).r_null_dist] = ...
            perm_rdm_corr(...
            res_adults.res(i_sub).rdm, ...
            models{model_id_newborns}, ...
            n_perm, ...
            'parfor', true, ...
            'covariate', rdm_stim);
    end
    save([fname_save_base, '_res_a2n_ind.mat'], 'res_a2n_ind');
end

%%
% RHO of individual subjects
f = figure('pos', [809 696 80 173], 'color', 'white');
pnl = panel(f);
ax = pnl.select();

plot_r_bar([res_a2n_ind.r_obs], [res_a2n_ind.p_perm], 'y_maxlim', 1, 'ax', ax);

ax.YLim = [-1, 1];
ax.YTick = [-1, 1];
ax.YAxis.LineWidth = 2;
ax.XAxis.LineWidth = 2;
txt = text(0, 1, sprintf('p=%.1g', res_a2n_group.p_perm));
txt.HorizontalAlignment = 'center';
pnl.margin = [8,4,0,1];

if fig_save
    saveas(f, [fname_save_base, '_res_adults_modelNewborns.png']);
end

% average RDMs across adults
for i_sub = 1:length(res_adults.res)
    rdm_adults_all(:, :, i_sub) = res_adults.res(i_sub).rdm;
    
end

f = plot_best_model(mean(rdm_adults_all, 3), models{model_id_newborns}, ...
    res_a2n_group.r_obs, res_a2n_group.p_perm, res_a2n_group.r_null_dist, ...
    'ratios', stim.ratios);
sgtitle('adults on newborn model');
if fig_save
    saveas(f, [fname_save_base, '_res_a2n.png']);
end


%% Hierarchical clustering
% plotting anchoring colours

datasets = {'newborns', 'adultsEEG'};

% re-define color depending on number of clusters: red, blue, orange,
% green, purple (with incremental inclusion
col_anchor = {{[0.8941, 0.1020, 0.1098], [0.2314, 0.4471, 0.6902]},
    {[0.8941, 0.1020, 0.1098], [1.0000, 0.4980, 0.0000], ...
    [0.2314, 0.4471, 0.6902]},
    {[0.8941, 0.1020, 0.1098], [1.0000, 0.4980, 0.0000], ...
    [0.3020, 0.6863, 0.2902], [0.2314, 0.4471, 0.6902]},
    {[0.8941, 0.1020, 0.1098], [1.0000, 0.4980, 0.0000], ...
    [0.5961, 0.3059, 0.6392], [0.3020, 0.6863, 0.2902], ...
    [0.2314, 0.4471, 0.6902]}};

n_cond = size(res_newborns.rdm.rdm, 1);

for i = 1:length(datasets)

    % compute dissimilarity matrix
    switch datasets{i}

        case 'newborns'

            DSM = 1 - res_newborns.rdm.rdm;

        case 'adultsEEG'

            % DSM from averaging RSMs across adult participants
            DSM = 1 - mean(rdm_adults_all, 3);
    end

    % make sure DSM is symmetric and contains 0 along the diagonal
    DSM = (DSM + DSM') / 2;
    DSM(1:size(DSM, 1)+1:end) = 0;

    % create hierarchical cluster tree
    % Z = linkage(DSM, 'average');
    Z = linkage(squareform(DSM), 'average');

    % define clusters
    clusters = [];
    clusters(:, 2) = cluster(Z, 'maxclust', 2);
    clusters(:, 3) = cluster(Z, 'maxclust', 3);
    clusters(:, 4) = cluster(Z, 'maxclust', 4);
    clusters(:, 5) = cluster(Z, 'maxclust', 5);

    % coordinates for 2-D multidimensional scaling
    [Y, stress] = cmdscale(DSM, 2);

    fhandle = figure('Position', [1 603 135 209]); % [404 555 218 314]
    f = panel();
    f.pack(size(clusters, 2)-1, 1);

    for i_cluster = 2 : size(clusters, 2)

        raw_labels = clusters(:, i_cluster); % raw cluster indices
        col_selection = col_anchor{i_cluster-1}; % ordered colors for this cluster

        % anchor cluster containing condition 1 to red and condition n_cond to
        % blue
        label_cond_start  = raw_labels(1); % raw label of condition 1
        label_cond_end = raw_labels(n_cond); % raw label of condition n_cond

        if isequal(label_cond_start, label_cond_end)
            warning(['same label for the two conditions at the extreme ' ...
                'of the continuum - plotting might not work properly']);
        end

        % create mapping between cluster label and color
        all_raw_labels = unique(raw_labels);
        col_map = zeros(1, i_cluster);
        col_map(label_cond_start) = 1;
        col_map(label_cond_end) = i_cluster;

        rest_col = setdiff(1:i_cluster, [1, i_cluster]);
        rest_labels = setdiff(all_raw_labels, [label_cond_start, label_cond_end]);

        for i_label = 1:length(rest_labels)
            col_map(rest_labels(i_label)) = rest_col(i_label);
        end

        f(i_cluster-1, 1).select();

        for i_cond = 1:n_cond
            col_idx = col_map(raw_labels(i_cond));
            text(Y(i_cond, 1), Y(i_cond, 2), num2str(i_cond), ...
                'color', col_selection{col_idx}, ...
                'FontWeight', 'bold', 'FontName', 'Arial', 'FontSize', 15);
            hold on;
        end

        xlim([min(Y(:, 1))-.2 max(Y(:, 1))+.2]);
        ylim([min(Y(:, 2))-.2 max(Y(:, 2))+.2]);
        set(gca, 'XTick', [], 'YTick', []);
        set(gca, 'LineWidth', 1.3);
    end

    f.title(sprintf('%s', datasets{i}));

    f.de.margin = 3;
    f.margin = [3 3 3 15];

    if fig_save
        print(fhandle, '-dpng', '-painters', '-r600', ...
            sprintf('%s_HC_%s_anchorcolour.png', fname_save_base, datasets{i}));
    end


    % plot dendrogram

    % specify dendrogram leaf order
    % D_for_leaf = pdist(DSM);
    % leaf_order = optimalleaforder(Z, D_for_leaf);
    leaf_order = optimalleaforder(Z, squareform(DSM));

    f = figure('Position', [476 653 266 217]);
    h_dend = dendrogram(Z, 'Reorder', leaf_order);
    set(h_dend, 'Color', 'k', 'LineWidth', 1.5);
    ax = gca;
    set(ax, 'FontName', 'Arial', 'Fontsize', 12);
    ylabel('dissimilarity', 'FontName', 'Arial', 'Fontsize', 15);
    xlabel('condition', 'FontName', 'Arial', 'Fontsize', 15);
    title(sprintf('%s', datasets{i}), 'FontName', 'Arial', 'Fontsize', 12);

    if fig_save
        print(f, '-dpng', '-painters', '-r600', ...
            sprintf('%s_dendrogram_%s.png', fname_save_base, datasets{i}));
    end
end
