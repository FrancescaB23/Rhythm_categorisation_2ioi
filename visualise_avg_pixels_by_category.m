%% Visualise average value per category

% script to plot the average values of the pixels belonging to different
% categories in newborns and adults (visualisation purposes)

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

mod_to_plot = 5; % number of models to plot

fig_save = true;


%%
fig_folder = fullfile(par_adults.deriv_path, 'figures', 'supplementary');


%%
harm_str = vec_to_str(harm_selections, 'format','%d', 'sep','-');

if feat_snr
    snr_str = sprintf('%d-%d', par_adults.snr_bins(1), par_adults.snr_bins(2));
else
    snr_str = 'none';
end

% save_path = fullfile(sprintf('%s_cat-2-3-4-5_comparisons', ...
%     par_adults.summary_path), 'avg_pixel_by_category_eeg');
% mkdir(save_path);


%% load results
% load results adults EEG
res_adults = load(fullfile(par_adults.summary_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf(['response-lf_feat-%s_harm-%s_snr-%s_' ...
    'roi-%s_avgChans-%s_cat-2-3-4-5_data.mat'], ...
    feat_resp, harm_str, snr_str, roi_name, jsonencode(avg_chans))));


% load results newborns
res_newborns = load(fullfile(par_newborns.results_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf('sub-grand-%s', ga_newborns), ...
    sprintf('sub-grand-%s_feat-%s_method-pearson_nCateg-2-3-4-5_perm.mat', ...
    ga_newborns, feat_resp)));

res_newborns.rdm = load(fullfile(par_newborns.results_path, ...
    sprintf('response-lf_stim-%s', feat_stim), ...
    sprintf('harm-%s', harm_str), ...
    sprintf('snr-%s', snr_str), ...
    sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)), ...
    sprintf('sub-grand-%s', ga_newborns), ...
    sprintf('sub-grand-%s_feat-%s_method-pearson_rdm.mat', ...
    ga_newborns, feat_resp)));


%% save fname

% fname_save_base = fullfile(save_path, ...
%     sprintf('response-lf_featStim-%s_harm-%s_concatEp-false_snr-%s_roi-%s_avgChans-%s_feat-%s_pixels-by-category_adultsVsNewborns', ...
%     feat_stim, ...
%     harm_str, ...
%     snr_str, ...
%     roi_name, ...
%     jsonencode(avg_chans), ...
%     feat_resp ...
%     ));


%% get model id sorted by highest correlation coefficient in newborns
% res_newborns.all_models_info

% sort result table by n of categories and categorical boundaries and
% assign model_id
res_newborns.all_models_info = sortrows(res_newborns.all_models_info, ...
    {'n_categ', 'bound1', 'bound2', 'bound3', 'bound4'});
res_newborns.all_models_info.mod_id = ...
    [1 : length(res_newborns.all_models_info.n_categ)]';

% get model id sorted by highest correlation coefficient (in newborns)
[~, model_id_descend] = sortrows(res_newborns.all_models_info, 'r', ...
    'descend');


%% categorical models
% generate all possible categorical models
n_categ = [2, 3, 4, 5];
n_mini_categ_skip = 1;
[models, all_models_info] = generate_categ_models(...
    par_adults.n_cond, n_categ, ...
    'n_mini_categ_skip', n_mini_categ_skip);


%% sanity check
% to make sure the model ordering is the same across the results in
% newborns and the generated categorical models

if (isequal(res_newborns.all_models_info.n_categ, all_models_info.n_categ) && ...
        isequal(res_newborns.all_models_info.bound1, all_models_info.bound1) && ...
        isequaln(res_newborns.all_models_info.bound2, all_models_info.bound2) && ...
        isequaln(res_newborns.all_models_info.bound3, all_models_info.bound3) && ...
        isequaln(res_newborns.all_models_info.bound4, all_models_info.bound4) && ...
        isequaln(res_newborns.all_models_info.bound5, all_models_info.bound5))
    disp('model order match');

else
    disp('WARNING model order does not match');
end


%% extract average pixel values for each category

% mask for values in the lower triangular part, excluding diagonal
mask = tril(ones(par_adults.n_cond, par_adults.n_cond, 'logical')) & ...
    ~diag(ones(1, par_adults.n_cond, 'logical')) ;
mask_idx = find(mask);


f = figure('Position', [1 314 1096 556]);


for i_mod = 1:mod_to_plot % models to plot


    % define pixels belonging to different categories
    cat_model = models{model_id_descend(i_mod)};

    % get a model matrix with numbers corresponding to category identity

    N = size(cat_model, 1);
    labels = zeros(N, 1);
    current_label = 0;


    for i = 1:N
        if labels(i) == 0  % not yet assigned
            current_label = current_label + 1;
            labels(i) = current_label;
        end
        % assign same label to all items connected to i
        connected = find(cat_model(i,:) == 1);
        for j = connected
            if labels(j) == 0
                labels(j) = labels(i);
            end
        end
    end


    % build the label matrix: value = category of row (or col)
    % within-category pairs get their category number, cross-category pairs get 0
    label_matrix = zeros(N, N);
    for i = 1:N
        for j = 1:N
            if i ~= j && labels(i) == labels(j)
                label_matrix(i,j) = labels(i);
            end
        end
    end

    % keep only values lower triangular part
    label_matrix(~mask) = 0;


    % extract values corresponding to each category and average them

    % get the number of categories
    categ_count = [];
    categ_count = unique(label_matrix);
    categ_count = categ_count(categ_count > 0);


    % loop over each category and extract average pixels
    avg_pixel_cat = [];
    for i_cat = 1:length(categ_count)

        mask_cat = label_matrix == i_cat; % category mask

        % newborns
        avg_pixel_cat.newborns(i_cat) = mean(res_newborns.rdm.rdm(mask_cat));

        % individual adults
        for i_sub = 1:length(res_adults.res)
            avg_pixel_cat.adults(i_sub, i_cat) = mean(res_adults.res(i_sub).rdm(mask_cat));
        end
    end


    % figure for supplementary material

    subplot(2, mod_to_plot, i_mod);
    
    bar(categ_count+.2, mean(avg_pixel_cat.adults), ...
        'FaceColor', [0 0 0], 'FaceAlpha', 0, ...
        'EdgeColor', [64 176 166]/255, 'LineWidth', 1.5);
    hold on;
    scatter(categ_count+.2, avg_pixel_cat.adults, 70, ...
        [64 176 166]/255, 'filled', 'MarkerFaceAlpha', 0.7, ...
        'MarkerEdgeColor', 'none');
    hold on;
    plot(categ_count-.2, avg_pixel_cat.newborns, '^', ...
        'MarkerFaceColor',[0 0 0], 'MarkerEdgeColor', 'none', ...
        'MarkerSize', 12);
    hold on;
    plot(categ_count-.2, avg_pixel_cat.newborns, '^', ...
        'MarkerFaceColor',[225 190 106]/255, 'MarkerEdgeColor', 'none', ...
        'MarkerSize', 10);
    yline(0, 'Color', 'k', 'LineWidth', 2);
    if max(categ_count) < 4;
        xlim([0.5-0.2 3.5+.2]);
    end
    ylim([-1.1 1.1]);

    set(gca, 'LineWidth', 2, 'YTick', [-1 0 1], 'XTick', 1:3, ...
        'TickLength', [0 0], 'Box', 'off', ...
        'FontName', 'Arial', 'FontSize', 13);
    ylabel('ρ', 'FontName', 'Arial', 'FontSize', 14);
    xlabel('category', 'FontName', 'Arial', 'FontSize', 14);
    title(sprintf('wm-%d id-%d ', i_mod, ...
        model_id_descend(i_mod)));

    subplot(2, mod_to_plot, i_mod + mod_to_plot);
    imagesc(cat_model.*labels);
    axis square;
    set(gca, 'XTick', [], 'YTick', [], 'TickLength', [0 0]);
    clim([0 3]); % additional, to keep color map consistent such that categories have the same colors regardless of the number of categories

    title(sprintf('wm-%d id-%d ', i_mod, ...
        model_id_descend(i_mod)));

end

if fig_save
    print(f, '-dpng', '-painters', '-r600', fullfile(fig_folder, ...
        sprintf('AVG-PIX_CAT_feat-%s_featStim-%s_harm-%s', ...
        feat_resp, feat_stim, harm_str)));
end