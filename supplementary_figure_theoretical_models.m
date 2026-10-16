clc
clear
close all


%% parameters
par = get_par('experiment', '2ioi_conga');

fig_save = true;


%%
figure_path = fullfile(par.deriv_path, 'figures', 'supplementary');
if ~exist(figure_path)
    mkdir(figure_path);
end


%% generate all possible categorical models (for plotting)
[models, all_models_info] = generate_categ_models_5(par.n_cond, [2, 3, 4, 5], ...
    'n_mini_categ_skip', par.n_mini_categ_skip);


%% plot all theoretical models (supplementary figure)
f = figure('Position', [1 533 656 333]); %[1 312 954 554]
t = tiledlayout(4, 9, 'TileSpacing', 'compact', 'Padding', 'compact');

for i = 1:length(models)
    ax = nexttile(t); % goes to row 2
    imagesc(ax, models{i});
    axis(ax,'off','square');
    title(sprintf('model %d', i), 'FontName', 'Arial', ...
        'FontWeight', 'bold', 'FontSize', 12);
end

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        fullfile(figure_path, 'theoretical_models_2-3-4-5-cat.png'));
end


%% correlate categorical models
% mask - lower triangular part of the matrix, diagonal excluded
mask = tril(ones(par.n_cond, par.n_cond, 'logical')) ...
    & ~diag(ones(1, par.n_cond, 'logical'));

for i = 1:length(models)

    for j = 1:length(models)

        model_i_lower_tri = models{i}(mask);
        model_j_lower_tri = models{j}(mask);

        [rho_across_models(i,j), p_across_models(i,j)] = corr(model_i_lower_tri, model_j_lower_tri, ...
            'Type', 'Spearman');

    end
end

f = figure('Position', [1 84 940 782]);
imagesc(rho_across_models);
axis equal tight;
% colorbar;

for i = 1:length(models)
    for j = 1:length(models)
        text(j, i, sprintf('%.2f', rho_across_models(i,j)), ...
            'HorizontalAlignment', 'center', ...
            'Color', 'black', 'FontName', 'Arial', ...
            'FontWeight', 'bold', 'FontSize', 8);
    end
end

labels = strcat("model ", string(1:33));
xticks(1:length(models));
xticklabels(labels);
yticks(1:length(models));
yticklabels(labels);
set(gca, 'TickLength', [0 0]);
set(gca, 'FontName', 'Arial');

if fig_save
    print(f, '-dpng', '-painters', '-r600', ...
        fullfile(figure_path, 'corr_across_theoretical_models_2-3-4-5-cat.png'));
end