function f = plot_r_over_all_models(r_all_models, tested_bounds, varargin)

parser = inputParser; 

addParameter(parser, 'ax', []); 

parse(parser, varargin{:}); 

ax = parser.Results.ax; 

f = []; 
if isempty(ax)
    f = figure('Color', 'white', 'Position', [1091 754 230 140]); 
    ax = axes(f); 
end

mu = mean(r_all_models, 1); 
sem = std(r_all_models, [], 1) ./ sqrt(size(r_all_models, 1)); 

plot(tested_bounds, mu, 'linew', 3)

hold on 
x = [tested_bounds, fliplr(tested_bounds)];
y = [mu + sem, fliplr(mu - sem)];
fill(x, y, 'b', 'FaceAlpha',0.1, 'EdgeColor', 'none');

box off
xlim([1, length(tested_bounds)]); 
xticks(tested_bounds); 
xlabel('bound')
ylabel('r')

