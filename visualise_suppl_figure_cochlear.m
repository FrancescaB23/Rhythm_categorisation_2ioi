%% Supplementary figure - cochlear model results
% Script to create the figures for the cochlear model results

clear
close all
clc


%% parameters
par = get_par('experiment', '2ioi_conga');

harmonics = [1:6];
feat_resp = 'mX';
feat_stim = 'ratio';

fig_save = true;


%%
fig_folder = fullfile(par.deriv_path, 'figures', 'supplementary');


%% load results - cochlear model responses in time
load_path = fullfile(par.deriv_path, 'cochlear_model');

an = cell(1, par.n_cond);

for i_cond=1:par.n_cond

    fname = sprintf('response-AN_cond-%02d.mat', i_cond);
    res = load(fullfile(load_path, fname));
    an{i_cond} = sum(res.AN.an_sout, 1);
end


%% plot responses over time - one pattern repetition - colour coded conditions
% define colours for the different conditions
cols = customcolormap([0, 1], {'#e87b0e', '#6a24a3'}, par.n_cond);

f = figure('color', 'white', 'position', [631 71 355 70*par.n_cond]);
pnl = panel(f);
pnl.pack('v', par.n_cond);

for i_cond=1:par.n_cond

    ax = pnl(i_cond).select();
    plot(an{i_cond}(1:floor(res.AN.Fs_an_sout*par.cycle_dur)), ...
        'color', cols(i_cond, :), ...
        'LineWidth', 2);
    axis(ax,'off');
    ax.YLim = [0, 500];

end

pnl.de.marginleft = 5;
pnl.de.marginright = 5;
pnl.de.margintop = 2;      % vertical space between tiles
pnl.de.marginbottom = 2;   % vertical space between tiles

if fig_save
    fname_save = 'COCHLEAR_urear-an_time-one-cycle_colour';
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_folder, fname_save), '.png']);
end


%% cochlear rdm
% load cochlear fft responses
fpath_cochlear = fullfile(['/Users/barbero/Documents/Projects/' ...
    'Rhototype_2ioi_newborns/adults/data/derivatives'], ...
    'fft', 'response-urear-an');
fname_cochlear = sprintf('response-urear-an_snr-none_%s', feat_resp);
[header_cochlear, data_cochlear] = CLW_load(fullfile(fpath_cochlear, fname_cochlear));

% get frequency vector - specific for cochlear responses
freq_cochlear = [0: header_cochlear.datasize(end)-1] * header_cochlear.xstep;
frex = 1 / par.cycle_dur * harmonics;

% get cochlear rdm
feat_for_rdm_cochlear = get_feat_from_fft(header_cochlear, data_cochlear, freq_cochlear, frex);

rdm_stim = get_rdm_from_feat(feat_for_rdm_cochlear, 'method', 'pearson');

f = figure('pos',[25 485 182 151], 'color', 'white');
imagesc(rdm_stim);
title(sprintf('cochlear %s', feat_resp));
axis off
axis image

if fig_save
    fname_save = sprintf('COCHLEAR_feat-%s_rdm', feat_resp);
    print(f, '-dpng', '-painters', '-r600', ...
        [fullfile(fig_folder, fname_save), '.png']);
end