function main_preproc_02_filter(par, sub)
% function to filter continuous data and segment it into trials
% newborn data: external ring of electrodes is removed
% adult data: keep only eeg channels


%%
% load
fpath = fullfile(par.deriv_path, 'merged', sub_num2str(sub)); 

fname = sprintf('sub-%03d_task-EEG', sub); 

[header, data] = CLW_load(fullfile(fpath, fname)); 


% arrange channels
switch par.experiment

    case 'newborn' % discard channels outer circle

        for ich = 1:length(par.chans_outer_circle)

            chans_outer_circle_idx(ich) = find( ...
                strcmp({header.chanlocs.labels}, ...
                par.chans_outer_circle(ich)));

        end

        data(:,chans_outer_circle_idx,:,:,:,:) = [];

        header.datasize = size(data);
        header.chanlocs(chans_outer_circle_idx) = [];

        
    case '2ioi_conga' % keep only eeg channels

        chan_labels = readtable('chan_labels.csv');
        [header, data] = RLW_arrange_channels(header, data, chan_labels.new);
end


% high-pass filter 
[header, data] = RLW_butterworth_filter(header, data, ...
                                        'filter_type', 'highpass', ...
                                        'low_cutoff', par.filt_cut_low, ...
                                        'filter_order', par.filt_order_low); 

% low-pass filter 
[header, data] = RLW_butterworth_filter(header, data, ...
                                        'filter_type', 'lowpass', ...
                                        'high_cutoff', par.filt_cut_high, ...
                                        'filter_order', par.filt_order_high); 


% segment all trials 
codes = unique({header.events.code}); 

segment_dur = par.trial_dur + 2 * par.epoch_buffer_dur; 

warning('segmenting trials with duration: %.2f s with %.2f s buffer on each side', ...
    par.trial_dur, par.epoch_buffer_dur); 

[header_ep, data_ep] = segment_safe(header, data, codes,...
                                        'x_start', -par.epoch_buffer_dur,...
                                        'x_duration', segment_dur); 


% save 
fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
mkdir(fpath); 

header_ep.name =  sprintf('sub-%03d_task-EEG_filt_ep', sub); 
CLW_save(fpath, header_ep, data_ep); 


end