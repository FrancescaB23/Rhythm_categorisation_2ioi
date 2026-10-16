% Preprocessing tapping data recorded from the tapping box

clc
clear
close all


%%
par = get_par('experiment', '2ioi_conga'); 

% subject number to process
sub = 20;


%%
tap_runs = [2, 4, 6]; 
    
datasets_tap = []; 
    
logs_tap = []; 
    
% prepare empty table where each row will be one tap trigger 
col_names = {'sub', 'run', 'trial', 'condition', 'trial_repeated', 'trial_duration', 'tap_onset'}; 
taps = cell2table(cell(0, length(col_names)), 'VariableNames', col_names); 

for run = tap_runs

    fprintf('run-%03d\n', run); 

    % load the logfile 
    d_log = dir(fullfile(par.raw_path, sub_num2str(sub), 'log', ...
                             sprintf('*run-%03d_*.tsv', run))); 

    if isempty(d_log)
        error('no log file found for %s', sub_num2str(sub)); 
    end
    
    log = readtable(fullfile(d_log.folder, d_log.name), ...
                        'Delimiter', '\t',...
                        'FileType', 'text', ...
                        'HeaderLines', 0, ...
                        'ReadVariableNames', true); 

    % remove bad trials (ones that were subsequently repeated)
    trial_events_to_keep_mask = ones(1, size(log, 1), 'logical'); 

    % check for any repeated trials (means the previous one was bad)
    for i_event=1:size(log, 1)
        if log{i_event, 'trial_repeated'}
            warning('removing bad trial event %d', i_event-1); 
            trial_events_to_keep_mask(i_event-1) = false; 
        end
    end

    % use the logfile to determine whether this is EEG or tapping block 
    is_tap_run = strcmp(log{1, 'trial_type'}{1}, 'TAP'); 

    % load source BDF file 
    d_bdf = dir(fullfile(par.raw_path, sub_num2str(sub), 'eeg',...
                         sprintf('*run-%03d_*.bdf', run))); 

    % load data 
    tic
    [header, data] =  RLW_import_BDF(fullfile(d_bdf.folder, d_bdf.name)); 
    fprintf('BDF import took %.1f seconds. \n', toc); 

    % only keep events with trigger value 1 or 2
    trial_events = header.events; 
    trial_events = trial_events(ismember({trial_events.code}, {'1', '2'})); 

    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % if there is an accidental, additional trigger, appearing for a
    % particular subject and run, you can set up a condition here to remove
    % it manually 
    if strcmp(par.experiment, '2ioi_conga') && sub == 5 && run == 6
        trial_events_to_keep_mask(10) = false;
        warning('sub 6 run 7 : keeping both trials for condition 1 in EEG file, event 10')
    end
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    % safety checks for trigger/log matching
    
    % make sure there is 
    if length(trial_events) ~= size(log, 1)
        error('there are %d events in the BDF file but %d in the logfile...', ...
              length(trial_events), size(log, 1))
    end

    % check that the inter-event intervals match (approx.) between BDF and
    % logfile 
    tap_trig_iois = diff([trial_events.latency])'; 
    log_file_trial_iois = diff(log.onset);

    disp([tap_trig_iois, log_file_trial_iois, tap_trig_iois-log_file_trial_iois]); 

    % this has to be a bit higher for 2ioi (5ms is fine for the other
    % experients)
    if any( abs(tap_trig_iois - log_file_trial_iois) > 0.020) 
        error('run-%03d: check trigger timing, there may be a problem!', run); 
    end

    % go over each event trigger 
    for i_event=1:length(trial_events)

        % rename the event code based on condition 
        new_code = sprintf('%d', log{i_event, 'condition_num'});
        trial_events(i_event).code = new_code; 

    end

    % ------------------------------------
    % keep the extracted data for this run 
    % ------------------------------------
    header.events = trial_events(trial_events_to_keep_mask);         

    if is_tap_run

        logs_tap = [logs_tap; log(trial_events_to_keep_mask, :)]; 

        datasets_tap = [datasets_tap, struct('header', header, 'data', data)]; 

    end
    
end


%% merge TAP datasets

[header, data] = merge_datasets(datasets_tap);

% make sure there are exactly 3 TAP runs
if length(datasets_tap) ~= 3
    error('there are only %d TAP runs instead of 3', length(datasets_tap));
end

% check one more time that we have exactly the same number of events in the
% lw file header and in the logflie for the merged dataset
if length(header.events) ~= size(logs_tap, 1)
    error('there are %d events in the BDF file but %d in the logfile...', ...
        length(header.events), size(logs_taps,1));
end


%% Prepare TAP datasets
% from continuous signal from the six force sensor to force signal

% segment all trials keeping 1 s before trial onset as baseline
codes = unique({header.events.code});

segment_dur = par.trial_dur + par.tap_baseline;

[header_ep, data_ep] = segment_safe(header, data, codes, ...
                            'x_start', -par.tap_baseline, ...
                            'x_duration', segment_dur);

% keep only force and trigger signal
channel_keep = ({header.chanlocs(80:86).labels});
[header_ep, data_ep] = RLW_arrange_channels(header_ep, data_ep, ...
                                        channel_keep);


%%
% make sure there are no taps in the baseline period
for trial = 1:size(data_ep, 1)

    f = figure('Position', [239 90 974 776]);

    for component = 1:6
        
        subplot(6, 1, component)
        plot(squeeze(data_ep(trial, component, 1, 1, 1, :)));
        xline(par.tap_baseline*1/header.xstep, 'r--', 'LineWidth', 3);
        hold on;
    end

    sgtitle(sprintf('trial %d', trial));
    drawnow;
    try 
        waitforbuttonpress
        close(f);
    catch
    end

end


%%
% baseline subtraction for each force component
header_f = header_ep;
data_f = data_ep;

for trial = 1:size(data_ep, 1)

    for component = 1:6
        
        % compute baseline 
        baseline = mean(data_ep(trial, component, 1, 1, 1, ...
            1:(par.tap_baseline * 1/header.xstep)));

        % baseline correction of force sensors
        components_corrected(component, :) = data_ep(trial, ...
                                component, 1, 1, 1, :) - baseline;

        % use calibration matrix to get from Volt to orthogonal force
        forces(component, :) = par.calib_coeff_Fz(component) .* components_corrected(component, :);
    end % end loop over components

    % Recompose force from components, pass from µV to V
    data_f(trial, size(data_ep,2)+1, 1, 1, 1, :) = -sum(forces)/1000000;

end % end loop over trials

header_f.datasize = size(data_f);
header_f.chanlocs(size(data_ep,2)+1).labels = 'Fz';
header_f.chanlocs(end).topo_enabled = 0;
header_f.chanlocs(end).SEEG_enabled = 0;

% segment to remove baseline
[header_f, data_f] = segment_safe(header_f, data_f, codes, ...
                            'x_start', 0, ...
                            'x_duration', par.trial_dur);

% downsample
[header_f, data_f] = RLW_downsample(header_f, data_f, ...
                                'x_downsample_ratio', 4);

% save file
fpath = fullfile(par.deriv_path, 'tapping', sub_num2str(sub));
mkdir(fpath);

header_f.name = sprintf('sub-%03d_task-TAP_response-force', sub);
CLW_save(fpath, header_f, data_f);