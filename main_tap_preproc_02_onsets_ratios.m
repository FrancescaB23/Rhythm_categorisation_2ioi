% Extract tap onset times from the continous data. This will depend on
% thresholds individually set for each subject so don't forget to stop the
% script manually and try different values, visually checking the
% extraction using the 'plot_diagnostic' parameter of extract_tap_ratios
% function.

clc
clear
close all


%% parameters
par = get_par('experiment', '2ioi_conga');

sub = 20


%% prepare output directories
sub_str = sub_num2str(sub);

save_path = fullfile(par.deriv_path, 'tap_onsets', sub_str);
mkdir(save_path);


%% load data
fpath = fullfile(par.deriv_path, 'tapping', sub_str);
fname = sprintf('sub-%03d_task-TAP_response-force.lw6', sub);
d = dir(fullfile(fpath, fname));

[header_ep, data_ep] = CLW_load(fullfile(d.folder, d.name));


%% load stimuli
load(fullfile(par.deriv_path, 'stimuli', 'stim.mat'));
stim_ratios = stim.ratios;
stim_iois = stim.ratios * stim.cycle_dur;


%%
fname = fullfile(save_path, 'onset_extraction_params.json');

if exist(fname)
    par_onset_extract = jsondecode(fileread(fname))
else
    par_onset_extract = [];
    par_onset_extract.tap_onset_thr = 0.16;% 0.10; % 0.07;
    par_onset_extract.tap_onset_min_asy = 0.030;
    par_onset_extract.tap_tone_min_asy = 0.124; %0.080;
    par_onset_extract.n_cycles_skip = 0;

    fid = fopen(fname, 'w');
    fprintf(fid, jsonencode(par_onset_extract));
    fclose(fid);
end

%%

tap_ratios_all = cell(stim.n_cond, 1);
n_valid_tap_ratios_all = zeros(stim.n_cond, 1);

for i_cond=1:par.n_cond

    % get only epochs for this condition
    ep_mask = strcmp({header_ep.events.code}, num2str(i_cond));
    [header_cond, data_cond] = RLW_arrange_epochs(header_ep, data_ep, find(ep_mask));

    % select only force channel

    [header_cond, data_cond] = RLW_arrange_channels( ...
        header_cond, data_cond, {'Fz'});


    for i_ep=1:header_cond.datasize(1)

        % extract tap times
        [tap_onset_times, tap_indices] = extract_taps(...
            ensure_row(squeeze(data_cond(i_ep, :, 1, 1, 1, :))), ...
            1/header_ep.xstep, ...
            par_onset_extract.tap_onset_thr, ...
            par_onset_extract.tap_onset_min_asy, ...
            'plot_diagnostic', true);

        % get ratios
        [tap_ratios, ~, n_valid_tap_ratios, n_possible_tap_ratios] = ...
            extract_tap_ratios(...
            tap_onset_times, ...
            stim_iois(i_cond, :), ...
            par.trial_dur / par.cycle_dur, ...
            'tap_tone_min_asy', par_onset_extract.tap_tone_min_asy, ...
            'n_cycles_skip', par_onset_extract.n_cycles_skip, ...
            'verbose', false, ...
            'plot_diagnostic', true);

        disp(length(tap_onset_times))

        drawnow;
        try
            waitforbuttonpress
            close(f);
        catch
        end

        close all % safety check if I accidentally leave plot_diagnostic=True

        tap_ratios_all{i_cond} = [tap_ratios_all{i_cond}; tap_ratios];
        n_valid_tap_ratios_all(i_cond) = n_valid_tap_ratios_all(i_cond) + n_valid_tap_ratios;

    end

end

tap_ratios_mean = cellfun(@(x) mean(x, 1), tap_ratios_all, 'uni', 0);
tap_ratios_mean = cat(1, tap_ratios_mean{:});


fname = sprintf('sub-%03d_response-tap-ratios.mat', sub);

tap_ratios = tap_ratios_mean;
n_valid_tap_ratios = n_valid_tap_ratios_all;
save(fullfile(save_path, fname), 'tap_ratios', 'stim_ratios', 'n_valid_tap_ratios');

n_valid_tap_ratios_all

rdm = get_rdm_from_feat(tap_ratios_mean', 'method', 'L2');

f = plot_rdm(rdm, ...
    'flip', false, ...
    'ratios', stim.ratios);
title(sprintf('sub-%03d',sub));

figure()
plot(stim_ratios(:,1), tap_ratios(:,1))
    