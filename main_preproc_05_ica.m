function main_preproc_05_ica(par, sub)
% function to run ICA


%%
% path to ICA matrix 
fpath_matrix = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
fname_matrix =  sprintf('sub-%03d_task-EEG_icaMatrix.mat', sub); 

% if the ICA matrix is saved, load it
if exist(fullfile(fpath_matrix, fname_matrix), 'file')
    
    load(fullfile(fpath_matrix, fname_matrix));
    warning('found ICA matrix on disk - loading that one!'); 
    
% if the ICA matrix is not saved, let's calcualte it
else

    %% load bads selected in the previous step 
    fname = fullfile(par.deriv_path, 'preproc', sub_num2str(sub), ...
                     sprintf('%s_badChans.csv', sub_num2str(sub))); 
    tbl_bad_chans = readtable(fname); 

    fname = fullfile(par.deriv_path, 'preproc', sub_num2str(sub), ...
                     sprintf('%s_badTrials.csv', sub_num2str(sub))); 
    tbl_bad_trials = readtable(fname); 

    bad_chans = tbl_bad_chans.channel;
    bad_epochs = tbl_bad_trials.trial; 

    %% preprocess data for ICA 
    fpath = fullfile(par.deriv_path, 'merged', sub_num2str(sub)); 
    fname = sprintf('sub-%03d_task-EEG', sub); 

    % load continuous data
    [header, data] = CLW_load(fullfile(fpath, fname)); 

    % high-pass filter
    [header, data] = RLW_butterworth_filter(header, data, ...
        'filter_type', 'highpass', ...
        'low_cutoff', par.filt_ica_cut_low, ...
        'filter_order', par.filt_ica_order_low);

    % low-pass filter
    [header, data] = RLW_butterworth_filter(header, data, ...
        'filter_type', 'lowpass', ...
        'high_cutoff', par.filt_ica_cut_high, ...
        'filter_order', par.filt_ica_order_high);

    % segment 
    codes = unique({header.events.code});
    [header_ep, data_ep] = segment_safe(header, data, codes,...
        'x_start', 0, 'x_duration', par.trial_dur);

    % downsample  
    [header_ep, data_ep] = RLW_downsample(header_ep, data_ep, ...
        'x_downsample_ratio', par.decim_factor_ica);

    % sanity check to prevent aliasing 
    if 1/header_ep.xstep < (4 * par.filt_ica_cut_high)
        warning('downsampling to %d Hz but antialiasing filter cutoff is only %d Hz',...
                1/header_ep.xstep,  par.filt_cut_high); 
    end


    % keep only EEG channels
    if strcmp(par.experiment, '2ioi_conga')
    
        chan_labels = readtable('chan_labels.csv'); % load channel labels
        [header_ep, data_ep] = RLW_arrange_channels(header_ep, data_ep, ...
            chan_labels.new);
    end


    %% reject bad trials 
    [header_ep, data_ep] = RLW_arrange_epochs(header_ep, data_ep, ...
        setdiff([1:header_ep.datasize(1)], bad_epochs));


    %% interpolate bad channels 
    if ~isempty(bad_chans)
        [header_ep, data_ep] = interpolate_bad_chans(header_ep, data_ep, bad_chans, ...
            par.n_chans_interp);
    end

    %% compute ICA 
    matrix_ICA = RLW_ICA_compute(header_ep, data_ep, ...
                                 'PICA_percentage', 99,...
                                 'ICA_mode', 'LAP');                           
    save(fullfile(fpath_matrix, fname_matrix), 'matrix_ICA');
    
end


%% ICA assign matrix
load('ica_gui_data.mat')

fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
fname =  sprintf('sub-%03d_task-EEG_icaMatrix', sub); 
load(fullfile(fpath, fname));

fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint', sub); 
header = CLW_load_header(fullfile(fpath, fname));

header.history = [];
header.history(end+1).configuration.gui_info = ica_gui_data;
header.history(end).configuration.parameters.ICA_um = matrix_ICA.ica_um;
header.history(end).configuration.parameters.ICA_mm = matrix_ICA.ica_mm;

% save the data back 
CLW_save_header(fpath, header);

end