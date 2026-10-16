function f = plot_boundary_distr(bounds, n_cond, varargin)

parser = inputParser; 

addParameter(parser, 'alpha', 0.1); 
addParameter(parser, 'pnl', []); 

parse(parser, varargin{:}); 

opacity = parser.Results.alpha; 
pnl = parser.Results.pnl; 

%%
f = []; 
if isempty(pnl)
    f = figure('color', 'white', 'Position',  [680 804 867 174]); 
    pnl = panel(f); 
end
pnl.pack('h', [30, 70]); 

rdms = cellfun(@(b) get_2categ_model(n_cond, b), num2cell(bounds), 'uni', 0); 

ax = pnl(1).select(); 

plot_rdms_overlay(rdms, 'ax', ax); 

%%

ax = pnl(2).select(); 

col = [0.6, 0.6, 0.6]; 
h = histogram(ax, bounds, [0.5 : 1 : n_cond-1+0.5], ...
    'EdgeColor', col, 'FaceColor', col, 'LineWidth', 2); 
box(ax, 'off'); 
ax.XTick = [1 : n_cond-1]; 
ax.YTick = [0 : 1 : ceil(ax.YLim(2))]; 
ax.FontSize = 12; 

