
clear all

par = get_par('experiment', '2ioi'); 

subjects = [6 : 23]; 

%%

save_path = fullfile(par.deriv_path, 'fft', 'sub-grand'); 

mkdir(save_path); 


%%

data_lf_all = cell(1, length(subjects)); 

for i_sub=1:length(subjects)

    sub = subjects(i_sub); 

    sub_str = sub_num2str(sub)

    % low-frequency response 
    % ----------------------

    % load data 
    fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
    fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt*.lw6', sub); 
    d = dir(fullfile(fpath, fname)); 

    [header_lf_ep, data_lf_ep] = CLW_load(fullfile(d.folder, d.name));

    % rereference 
    [header_lf_ep, data_lf_ep] = RLW_rereference(header_lf_ep, data_lf_ep, ...
                                    'apply_list', {header_lf_ep.chanlocs.labels}, ...
                                    'reference_list', par.eeg_ref_chans); 

    % don't concatenate epochs
    [header_lf_avg, data_lf_avg] = average_trials_per_cond(header_lf_ep, data_lf_ep, ...
                                'concatenate_epochs', false); 
                    
    data_lf_all{i_sub} = data_lf_avg; 
    

end



%%

header_lf = header_lf_avg; 

data_lf = mean(cat(7, data_lf_all{:}), 7); 

fname = sprintf('sub-grand_response-lf_concatEp-false'); 

get_and_save_ffts(header_lf, data_lf, par.snr_bins, ...
                  'fname', fname, 'save_path', save_path); 

              

