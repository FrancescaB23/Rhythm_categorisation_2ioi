function main_tap_preproc_05_erp(par, sub)
% function to compute erp of tap data, all chunks included


%% prepare output directory
sub_str = sub_num2str(sub); 
save_path = fullfile(par.deriv_path, 'erp', sub_str); 

if ~exist(save_path)
    mkdir(save_path);
end


%% load and preprocess data
% load data
fpath = fullfile(par.deriv_path, 'tapping', sub_num2str(sub));
fname =  sprintf('sub-%03d_task-TAP_response-force', sub);

[header, data] = CLW_load(fullfile(fpath, fname));

% keep only Fz channel
[header, data] = RLW_arrange_channels(header, data, {'Fz'});


%%
% low-pass filter for ERPs
[header, data] = RLW_butterworth_filter(header, data, ...
    'filter_type', 'lowpass', ...
    'high_cutoff', par.filt_tap_erp_cut_low, ...
    'filter_order', 4);


% re-segment
codes = unique({header.events.code});
[header, data] = segment_safe(header, data, codes,...
    'x_start', 0, 'x_duration', par.trial_dur);


% get erp chunks
[header_chunk, data_chunk] = chunk_cycles_per_cond(...
    header, data, par.cycle_dur, par.experiment);


% downsample data to 64 Hz 
[header_chunk, data_chunk] = RLW_downsample(header_chunk, data_chunk, ...
    'x_downsample_ratio', par.decim_factor);
if (1/header_chunk.xstep) ~= 64
    error('error: check sampling rate and adjust downsampling accordingly');
end

                           
% baseline correction (subtract mean of each epoch)
data_chunk = data_chunk - mean(data_chunk, 6); 


% average cycle chunks separately for each condition
[header_erp, data_erp] = average_trials_per_cond( ...
    header_chunk, data_chunk, par, 'concatenate_epochs', false);
     

% save erp data
fname = sprintf('sub-%03d_response-tap-force-erp_lpf-%.0fHz', ...
    sub, par.filt_erp_cut_low); 
    
header_erp.name = fname; 
CLW_save(save_path, header_erp, data_erp); 


% sanity check plot 
plot_erp_overlay(squeeze(mean(data_erp, 2)), 1/header_erp.xstep);  
title(sub_str);
end