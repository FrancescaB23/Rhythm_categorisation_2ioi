
clear

subjects = [6 : 23]; 

par = get_par('experiment', '2ioi'); 

for i_sub=1:length(subjects)
    
    sub = subjects(i_sub); 
    sub_str = sub_num2str(sub); 
    
    % prepare output directories
    save_path = fullfile(par.deriv_path, 'fft', sub_str); 
    mkdir(save_path); 

    % load     
    fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
    fname = sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt', sub); 

    [header_ep, data_ep] = CLW_load(fullfile(fpath, fname)); 
        
    % rereference 
    [header_ep, data_ep] = RLW_rereference(header_ep, data_ep, ...
                                    'apply_list', {header_ep.chanlocs.labels}, ...
                                    'reference_list', par.eeg_ref_chans); 

    % don't concatenate epochs
    [header_avg, data_avg] = average_trials_per_cond(header_ep, data_ep, ...
                                'concatenate_epochs', false); 

    fname = sprintf('sub-%03d_response-lf_concatEp-false', sub); 

    get_and_save_ffts(header_avg, data_avg, par.snr_bins, ...
                      'fname', fname, 'save_path', save_path, ...
                      'fmax', 64); 

    % concatenate epochs
    [header_avg, data_avg] = average_trials_per_cond(header_ep, data_ep, ...
                                'concatenate_epochs', true); 

    fname = sprintf('sub-%03d_response-lf_concatEp-true', sub); 

    get_and_save_ffts(header_avg, data_avg, par.snr_bins, ...
                      'fname', fname, 'save_path', save_path, ...
                      'fmax', 64); 


end
