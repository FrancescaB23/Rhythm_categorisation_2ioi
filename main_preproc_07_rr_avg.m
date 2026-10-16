function main_preproc_07_rr_avg(par, sub)
% function to further process data after artefact rejection:
% - remove some more external electrodes for consistency with Edalati et
% al. 2023 (newborn data only)
% - segment [0; trial duration]
% - downsample
% - rereference
% - average trials per condition


%%
% load data
fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub));
if strcmp(par.experiment, 'newborn')
    fname =  sprintf('sub-%03d_task-EEG_filt_ep_AB-100_chanint.lw6', sub);
else
    fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*.lw6', sub);
end
d = dir(fullfile(fpath, fname));

[header, data] = CLW_load(fullfile(d.folder, d.name));


% newborns data only: remove some more 'external' channels for consistency 
% with Edalati et al. 2023
if strcmp(par.experiment, 'newborn')

    chans_outer_ring = {'E14' 'E17' 'E21' 'E43' 'E48' 'E38' 'E44' ...
        'E49' 'E56' 'E107' 'E113' 'E120' 'E1' 'E32' 'E25' 'E8' ...
        'E121' 'E114' 'E99' 'E94' 'E88' 'E81' 'E73' 'E68' 'E63' ...
        'E119' 'E125' 'E126' 'E127' 'E128'};

    chans_outer_ring_idx = find(ismember( ...
        {header.chanlocs.labels}, chans_outer_ring));

    data(:,chans_outer_ring_idx,:,:,:,:) = [];

    header.datasize = size(data);
    header.chanlocs(chans_outer_ring_idx) = [];
end


% re-segment
codes = unique({header.events.code});
[header, data] = segment_safe(header, data, codes,...
    'x_start', 0, 'x_duration', par.trial_dur);


% downsample
[header, data] = RLW_downsample(header, data, ...
    'x_downsample_ratio', par.decim_factor);


% rereference to average mastoids
[header, data] = RLW_rereference(header, data, ...
    'apply_list', {header.chanlocs.labels}, ...
    'reference_list', par.eeg_ref_chans);


% average trials per condition
[header_avg, data_avg] = average_trials_per_cond(header, data, par, ...
    'concatenate_epochs', false);

header_avg.name = [header_avg.name, '_rr-m'];


% save files
CLW_save(fpath, header_avg, data_avg);

end