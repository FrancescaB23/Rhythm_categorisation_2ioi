function [f,pnl] = plot_time_mX_sim(x_all, fs, frex, varargin)

parser = inputParser(); 

addParameter(parser, 'pnl', []); 
addParameter(parser, 'cycle_dur', 0.750); 
addParameter(parser, 'filter_erp', true); 
addParameter(parser, 'cond_labs', []); 

parse(parser, varargin{:}); 

pnl = parser.Results.pnl; 
cycle_dur = parser.Results.cycle_dur; 
filter_erp = parser.Results.filter_erp; 
cond_labs = parser.Results.cond_labs; 

%%

if iscell(x_all)
    x_all = cat(1, x_all{:}); 
end

N = size(x_all, 2); 

n_cond = size(x_all, 1); 

[X, freq] = get_X(x_all, fs); 

mX = abs(X) ./ N * 2; 


%%

f = []; 
if isempty(pnl)
    f = figure('color', 'white', 'position', [631 71 600 930]); 
    pnl = panel(f); 
end

pnl.pack('h', [30, 70]); 
pnl(1).pack('v', repmat({[]}, 1, n_cond)); 
pnl(2).pack('v', repmat({[]}, 1, n_cond)); 

pnl.de.margin = 0;

fmin = 0.1; 
fmax = max(max(frex), 15); 
y_min = Inf; 
y_max = -Inf;


for i_cond=1:n_cond
    
    if filter_erp
        [b, a] = butter(2, 10/(fs/2)); 
        x = filtfilt(b, a, x_all(i_cond, :)); 
    else
        x = x_all(i_cond, :); 
    end
    x = epoch_chunks(x, fs, cycle_dur); 
    x = mean(x, 1); 
    t = [0 : length(x)-1] / fs; 

    frex_idx = round(frex / fs * N) + 1; 
    min_freq_idx = round(fmin / fs * N) + 1; 
    max_freq_idx = round(fmax / fs * N) + 1; 

    mX_to_plot = mX(i_cond, :); 
    mX_to_plot(1) = 0; 

    ax = pnl(1, i_cond).select(); 

    plot(t, x, 'k', 'linew', 1.5); 
    hold(ax, 'on'); 

    ax.XLim = [0, cycle_dur * 1]; 
    ax.XTick = [0, cycle_dur * 1]; 
    ax.YTick = [];     

    if i_cond ~= n_cond
        ax.XTick = []; 
        ax.YTick = []; 
    end

    ax = pnl(2, i_cond).select(); 

    plot(freq(min_freq_idx:max_freq_idx), ...
         mX_to_plot(min_freq_idx:max_freq_idx), ...
         'linew', 1); 

    hold on 

    plot(freq(frex_idx), mX_to_plot(frex_idx), 'ro'); 

    if i_cond ~= n_cond
        ax.XTick = []; 
        ax.YTick = []; 
        ax.XAxis.Visible = 'off'; 
    end

    y_min = min(y_min, min(mX_to_plot(min_freq_idx:max_freq_idx))); 
    y_max = max(y_max, max(mX_to_plot(min_freq_idx:max_freq_idx))); 

    if ~isempty(cond_labs)
        tit = pnl(1, i_cond).title(cond_labs{i_cond}); 
    end
    
end

for i_cond=1:n_cond
    ax = pnl(2, i_cond).select(); 
    ax.YLim = [y_min, y_max]; 
    ax.Title.Position(2) = ax.Title.Position(2) - ax.Title.Position(2) * 0.3; 
end

ax.YTick = [y_min, y_max]; 
ax.YTickLabel = [round(y_min * 100)/100, round(y_max * 100)/100]; 
ax.TickDir = 'out'; 

pnl.de.margintop = 8;
pnl.de.marginleft = 15; 
pnl.marginbottom = 10; 
pnl.margintop = 10; 
pnl.fontsize = 12; 
