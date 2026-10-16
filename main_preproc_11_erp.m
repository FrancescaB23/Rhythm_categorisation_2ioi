function main_preproc_11_erp(par, sub)
% function to compute erp of eeg data, all chunks included


%% prepare output directory
sub_str = sub_num2str(sub); 
save_path = fullfile(par.deriv_path, 'erp', sub_str); 

if ~exist(save_path)
    mkdir(save_path);
end


%% load and preprocess data
% load data
fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub));
fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*.lw6', sub);
d = dir(fullfile(fpath, fname));

% discard file with rereferenced data (as they are already averaged across
% trials per condition)
d = d(~contains({d.name}, 'rr-m', 'IgnoreCase', true));

% make sure just one file is selected, otherwise error
if numel(d) > 1
    error('sub-%03d: expected 1 file, found %d:\n%s', sub, numel(d), strjoin({d.name}, newline));
end

d.name
[header, data] = CLW_load(fullfile(d.folder, d.name));


%%
% low-pass filter for ERPs
[header, data] = RLW_butterworth_filter(header, data, ...
    'filter_type', 'lowpass', ...
    'high_cutoff', par.filt_erp_cut_low, ...
    'filter_order', 4);


% re-segment
codes = unique({header.events.code});
[header, data] = segment_safe(header, data, codes,...
    'x_start', 0, 'x_duration', par.trial_dur);


% rereference to average mastoids
[header, data] = RLW_rereference(header, data, ...
    'apply_list', {header.chanlocs.labels}, ...
    'reference_list', par.eeg_ref_chans);


% get erp chunks
[header_chunk, data_chunk] = chunk_cycles_per_cond(...
    header, data, par.cycle_dur, par.experiment);


% downsample data to 64 Hz (factor 16)
[header_chunk, data_chunk] = RLW_downsample(header_chunk, data_chunk, ...
    'x_downsample_ratio', 16);
if (1/header_chunk.xstep) ~= 64
    error('error: check sampling rate and adjust downsampling accordingly');
end

                           
% baseline correction (subtract mean of each epoch)
data_chunk = data_chunk - mean(data_chunk, 6); 


% average cycle chunks separately for each condition
[header_erp, data_erp] = average_trials_per_cond( ...
    header_chunk, data_chunk, par, 'concatenate_epochs', false);
     

% save erp data
fname = sprintf('sub-%03d_response-erp_lpf-%.0fHz_rr-m', ...
    sub, par.filt_erp_cut_low); 
    
header_erp.name = fname; 
CLW_save(save_path, header_erp, data_erp); 


% sanity check plot 
plot_erp_overlay(squeeze(mean(data_erp, 2)), 1/header_erp.xstep);  
title(sub_str);
end