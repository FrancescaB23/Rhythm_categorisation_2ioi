clear 

par = get_par('experiment', '2ioi'); 

load('/Users/tomaslenc/projects_backed_up/Rhytotype/data_eeg_2ioi/derivatives/tapping_francesca/TAP-force_18-sub.mat'); 

%%

load(fullfile(par.deriv_path, 'stimuli', 'stim.mat')); 

stim_ratios = stim.ratios; 
stim_iois = stim.ratios * stim.cycle_dur; 


%% 

fs = 1024; 
xstart = -1; 

for i_sub=1:length(TAP)
    
    tap_data = TAP{i_sub}; 
    
    sub = tap_data.subject_code; 
    sub_str = sub_num2str(sub); 
    
    sub_str
    
    save_path = fullfile(par.deriv_path, 'tap_onsets', sub_str); 
    mkdir(save_path); 

    fname = fullfile(save_path, 'onset_extraction_params.json');    
    if exist(fname)
        par_onset_extract = jsondecode(fileread(fname)); 
    else
        par_onset_extract = []; 
        par_onset_extract.tap_onset_thr = 0.1; 
        par_onset_extract.tap_onset_min_asy = 0.030; 
        par_onset_extract.tap_tone_min_asy = 0.080; 
        par_onset_extract.n_cycles_skip = 0; 

        fid = fopen(fname, 'w'); 
        fprintf(fid, jsonencode(par_onset_extract)); 
        fclose(fid); 
    end
       
    data = nan(13*3, 1, 1, 1, 1, round(par.trial_dur * fs)); 
    events = []; 
    
    tap_ratios_all = cell(stim.n_cond, 1); 
    n_valid_tap_ratios_all = zeros(stim.n_cond, 1); 
    
    c = 1; 
    for i_cond=1:13        
        for i_trial=1:3
            
            idx_start = round( (0 - xstart) * fs) + 1; 
                        
            trigger_on = ensure_row(tap_data.data{i_cond, i_trial}.trigger_on); 
            trigger_on = trigger_on - min(trigger_on); 
            trigger_on = trigger_on ./ max(trigger_on); 
            trigger_on = trigger_on(idx_start : end); 
            
            fz = ensure_row(tap_data.data{i_cond, i_trial}.Fz); 
            fz = fz(idx_start : end); 

            s = ensure_row(tap_data.data{i_cond, i_trial}.stimuli); 
            s = s - min(s); 
            s = s ./ max(s); 
            s = s(idx_start : end); 
                       
            % get tap onset times
            [tap_onsets, tap_index] = extract_taps(...
                                           trigger_on, ...
                                           fs, ...
                                           par_onset_extract.tap_onset_min_asy, ...
                                           par_onset_extract.tap_onset_min_asy, ...
                                           'plot_diagnostic', false); 

            % get ratios 
            [tap_ratios, ~, n_valid_tap_ratios, n_possible_tap_ratios] = ...
                        extract_tap_ratios(...
                                tap_onsets, ...
                                stim_iois(i_cond, :), ...
                                par.trial_dur / par.cycle_dur, ...
                                'tap_tone_min_asy', par_onset_extract.tap_tone_min_asy, ...
                                'n_cycles_skip', par_onset_extract.n_cycles_skip, ...
                                'plot_diagnostic', false); 

%             %t = [0 : length(trigger_on)-1] / fs - xstart; 
%             t = [0 : length(trigger_on)-1] / fs; 
% 
%             figure('Position', [620 727 1162 320])
%             plot(t, trigger_on, 'linew', 2); 
%             hold on 
%             plot(t, s); 
%             fz_norm = fz; 
%             fz_norm = fz_norm ./ max(fz_norm); 
%             fz_norm = fz_norm - min(fz_norm); 
%             plot(t, fz_norm, 'color', [.5 .5 .5]); 
%             plot([tap_onsets', tap_onsets'], [-1, 2], 'm--')
%             xlim([-1, 23]); 
                                   
            data(c, 1, 1, 1, 1, :) = fz; 
            
            events(c).code = num2str(i_cond); 
            events(c).latency = 0; 
            events(c).epoch = c; 

            tap_ratios_all{i_cond} = [tap_ratios_all{i_cond}; tap_ratios]; 
            n_valid_tap_ratios_all(i_cond) = n_valid_tap_ratios_all(i_cond) + n_valid_tap_ratios; 
                        
            c = c+1; 
            
        end        
    end
    
    % get mean tap ratio
    tap_ratios_mean = cellfun(@(x) mean(x, 1), tap_ratios_all, 'uni', 0); 
    tap_ratios_mean = cat(1, tap_ratios_mean{:}); 

    % save tap ratios 
    fname = sprintf('sub-%03d_response-tap-ratios.mat', sub); 
    tap_ratios = tap_ratios_mean; 
    n_valid_tap_ratios = n_valid_tap_ratios_all; 
    save(fullfile(save_path, fname), 'tap_ratios', 'stim_ratios', 'n_valid_tap_ratios'); 

    % plot RDM (sanity check)
    rdm = get_rdm_from_feat(tap_ratios_mean', 'method', 'L2'); 
    f = plot_rdm(rdm, ...
                 'flip', true, ...
                 'ratios', stim.ratios); 

    % prepare header
    header = make_lw_header(...
                        'data', data, ...
                        'fs', fs, ...
                        'events', events); 

    header.chanlocs.labels = 'force'; 
    
    % save continuous data in lw6 format
    fpath = fullfile(par.deriv_path, 'tapping', sub_str); 
    mkdir(fpath); 
    header.name = sprintf('sub-%03d_task-TAP_response-force', sub);
    CLW_save(fpath, header, data); 
    
end

