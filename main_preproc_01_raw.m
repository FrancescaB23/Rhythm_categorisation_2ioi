% This script does kind of source2raw conversion 
% Here, we do some cleaning and moving of files, we go over the runs
% and make sure all the triggers are okay, concatenate all the runs into a
% single continuous file and save it to derivatives. 

clear 
clc
par = get_par('experiment', '2ioi_conga'); % 'newborn' '2ioi_conga'

% subject number to process
sub = 20;  

%%

if strcmp(par.experiment, 'newborn') && sub ~= 2

    % load the logfile
    d_log = dir(fullfile(par.raw_path,sprintf('M%02d',sub),'log*.tsv'));
    
    if isempty(d_log)
        error('no log file found for %s', sub_num2str(sub)); 
    end
    
    log = readtable(fullfile(d_log.folder, d_log.name), ...
                    'Delimiter', '\t',...
                    'FileType', 'text', ...
                    'HeaderLines', 0, ...
                    'ReadVariableNames', true); 
    
    % load source file
    d_mff = fullfile(par.raw_path,sprintf('M%02d',sub), ...
        sprintf('M%02d.mff',sub));
    
    % load data
    tic
    [header, data] = RLW_import_EGI(d_mff);
    fprintf('MFF import took %.1f seconds. \n', toc); 
    
    
    %% Trigger checks:
    % 1- check match between triggers in log and eeg
    % 2- check inter-events intervals match between log and eeg (approximately)
    % 3- delete repeated trials
    
    % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    % if there is an accidental, additional trigger, appearing for a
    % particular subject and run, you can set up a condition here to remove
    % it manually 
    if sub == 1

        header.events(end) = []; 
        log(1,:) = [];
        warning('sub 1 removing last trial from EEG file and first trial from log file as corresponding trigger is missing in EEG');

    elseif sub == 3

        header.events(end) = []; 
        log(1,:) = [];
        % forcing this value to 0 as, even if the trial has been repeated, 
        % the trigger from the previous trial is not present in the EEG file
        log.trial_repeated(1) = 0; 
        warning(['sub 3 removing last trial from EEG file as it has a very long latency and it is not present in the logfile. ' ...
            'Remove first trial from log file as corresponding trigger is missing in EEG']);

    end
    
    % Remove from the header triggers not related to trial onset
    header.events = header.events(startsWith({header.events.code},'S'));
    
    % Retrieve triggers' code and latency 
    eeg_trigger = header.events;
    
    for i = 1:length(eeg_trigger)
        eeg_trigger_code(i) = str2num(eeg_trigger(i).code(end-1:end));
        eeg_trigger_latency(i) = eeg_trigger(i).latency;
    end
    
    % Check match between triggers from log and eeg 
    disp([eeg_trigger_code', log.condition_num])
    if ~isequal(eeg_trigger_code',log.condition_num)
        error(['Mismatch trigger codes between eeg and log file, there may be a problem']);
    end
    
    % Check that the inter-event intervals match (approx.) between EEG and log
    eeg_trigger_ioi = diff(eeg_trigger_latency)';
    log_trigger_ioi = diff(log.onset);
    
    disp([eeg_trigger_ioi, log_trigger_ioi, eeg_trigger_ioi - log_trigger_ioi])
    
    if any(abs(eeg_trigger_ioi - log_trigger_ioi) > 0.021)
        error('Check trigger timing, there may be a problem!');
    end
    
    % Delete repeated trials and aborted trials
    trial_events_to_keep_mask = ones(1, size(log, 1), 'logical'); 
    
    
    for i_event=1:size(log, 1)

        if log{i_event, 'trial_repeated'}
            warning('removing bad trial event %d', i_event-1); 
            trial_events_to_keep_mask(i_event-1) = false; 
        end
        
        if log{i_event, 'trial_aborted'}
            warning('removing bad trial event %d', i_event);
            trial_events_to_keep_mask(i_event) = false;
        end
    
    end
    
    header.events = header.events(trial_events_to_keep_mask);



elseif strcmp(par.experiment, 'newborn') && sub == 2

    % load the logfile
    d_log = dir(fullfile(par.raw_path,sprintf('M%02d',sub),'log*.tsv'));
    
    if isempty(d_log)
        error('no log file found for %s', sub_num2str(sub)); 
    end
    
    log = readtable(fullfile(d_log.folder, d_log.name), ...
                    'Delimiter', '\t',...
                    'FileType', 'text', ...
                    'HeaderLines', 0, ...
                    'ReadVariableNames', true); 

    % Remove first trial from the log file as the corresponding trigger is
    % missing in EEG
    log(1,:) = [];
    warning('sub 2 removing first trial from log file as corresponding trigger is missing in EEG');

    datasets_eeg = []; 


    for ifile = [1 2] % data for M02 are split into two files
    
        % load source file
        d_mff = fullfile(par.raw_path,sprintf('M%02d',sub), ...
            sprintf('M%02d_%d.mff',sub,ifile));
    
        % load data
        tic
        [header, data] = RLW_import_EGI(d_mff);
        fprintf('MFF import took %.1f seconds. \n', toc); 
    
    
        %% Trigger checks:
        % 1- check match between triggers in log and eeg
        % 2- check inter-events intervals match between log and eeg (approximately)
        % 3- delete repeated trials
        
        % !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        % if there is an accidental, additional trigger, appearing for a
        % particular subject and run, you can set up a condition here to remove
        % it manually 

        % Remove from the header triggers not related to trial onset
        header.events = header.events(startsWith({header.events.code},'S'));
        log_idx = [];

        if ifile == 1
            header.events(end-1:end) = []; 
            warning('sub 2, file 1: removing last two triggers from EEG file as trials were aborted');

            % Indices for correspondance between eeg triggers and trials in the
            % log file (as there is one logfile for two eeg files
            log_idx = [1:length(header.events)];

        elseif ifile == 2

            % Indices for correspondance between eeg triggers and trials in the
            % log file (as there is one logfile for two eeg files
            log_idx = [size(log,1)-length(header.events)+1:size(log,1)];

        end
    

        % Retrieve triggers' code and latency 
        eeg_trigger = header.events;
        
        for i = 1:length(eeg_trigger)
            eeg_trigger_code(i) = str2num(eeg_trigger(i).code(end-1:end));
            eeg_trigger_latency(i) = eeg_trigger(i).latency;
        end


        % Check match between triggers from log and eeg 
        if ~isequal(eeg_trigger_code',log.condition_num(log_idx))
            error(['Mismatch trigger codes between eeg and log file, thre may be a problem']);
        end
    
        % Check that the inter-event intervals match (approx.) between EEG and log
        eeg_trigger_ioi = diff(eeg_trigger_latency)';
        log_trigger_ioi = diff(log.onset(log_idx));
    
        disp([eeg_trigger_ioi, log_trigger_ioi])
    
        if any(abs(eeg_trigger_ioi - log_trigger_ioi) > 0.020)
            error('Check trigger timing, there may be a problem!');
        end
    
        % Delete repeated trials and aborted trials - done manually for M02


        datasets_eeg = [datasets_eeg, struct('header', header, 'data', data)]; 

    end

    % merge EEG datasets
    
    [header, data] = merge_datasets(datasets_eeg);


    
else

    eeg_runs = [1, 3, 5, 7]; 
    
    datasets_eeg = []; 
    
    logs_eeg = []; 
    
    % prepare empty table where each row will be one tap trigger 
    col_names = {'sub', 'run', 'trial', 'condition', 'trial_repeated', 'trial_duration', 'tap_onset'}; 
    taps = cell2table(cell(0, length(col_names)), 'VariableNames', col_names);  
    
    for run = eeg_runs
    
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
        is_eeg_run = strcmp(log{1, 'trial_type'}{1}, 'EEG'); 
    
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
        if strcmp(par.experiment, '2ioi_conga') && sub == 6 && run == 1
            trial_events_to_keep_mask(1) = false;
            warning('sub 6 run 1 : removing trial 1 from EEG file')
        end
        if strcmp(par.experiment, '2ioi_conga') && sub == 6 && run == 7
            trial_events_to_keep_mask(10) = true;
            warning('sub 6 run 7 : keeping both trials for condition 1 in EEG file, event 10')
        end
        if strcmp(par.experiment, '2ioi_conga') && sub == 17 && run == 3
            trial_events_to_keep_mask(11) = false;
            trial_events_to_keep_mask(12) = false;
            trial_events_to_keep_mask(13) = false;
            trial_events_to_keep_mask(14) = false;
            trial_events_to_keep_mask(15) = false;
            warning('sub 17 run 3 : removing trial 11 from EEG file')
            warning('sub 17 run 3 : removing trial 12 from EEG file')
            warning('sub 17 run 3 : removing trial 13 from EEG file')
            warning('sub 17 run 3 : removing trial 14 from EEG file')
            warning('sub 17 run 3 : removing trial 15 from EEG file')
        end
        if strcmp(par.experiment, '2ioi_conga') && sub == 19 && run == 5
            trial_events_to_keep_mask(10) = false;
            trial_events_to_keep_mask(11) = false;
            trial_events_to_keep_mask(12) = false;
            trial_events_to_keep_mask(13) = false;
            warning('sub 19 run 5 : removing trial 10 from EEG file')
            warning('sub 19 run 5 : removing trial 11 from EEG file')
            warning('sub 19 run 5 : removing trial 12 from EEG file')
            warning('sub 19 run 5 : removing trial 13 from EEG file')
        end

        if strcmp(par.experiment, '2ioi_conga') && sub == 19 && run == 7
            trial_events_to_keep_mask(:) = true;
            warning('sub 19 run 7 : keeping all trials')
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
        eeg_trig_iois = diff([trial_events.latency])'; 
        log_file_trial_iois = diff(log.onset);

        % this match cannot be done for trials occurring after a pause in
        % the recordings - check manually
        if strcmp(par.experiment, '2ioi_conga') && sub == 17 && run == 3
            eeg_trig_iois(9) = [];
            log_file_trial_iois(9) = [];
        end
    
        disp([eeg_trig_iois, log_file_trial_iois, eeg_trig_iois - log_file_trial_iois]); 
    
        % this has to be a bit higher for 2ioi (5ms is fine for the other
        % experients)
        if any( abs(eeg_trig_iois - log_file_trial_iois) > 0.020) 
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
    
        if is_eeg_run
    
            logs_eeg = [logs_eeg; log(trial_events_to_keep_mask, :)]; 
    
            datasets_eeg = [datasets_eeg, struct('header', header, 'data', data)]; 
    
        end
    
    end
    
    
    %% merge EEG datasets
    
    [header, data] = merge_datasets(datasets_eeg);
    
    % make sure there are exactly 4 EEG runs  
    if length(datasets_eeg) ~= 4
        error('there are only %d EEG runs instead of 4', length(datasets_eeg)); 
    end
    
    % check one more time that we have exactly the same number of events in the lw
    % file header and in the logfile (for the merged dataset)
    if length(header.events) ~= size(logs_eeg, 1)
        error('there are %d events in the BDF file but %d in the logfile...', ...
              length(header.events), size(logs_eeg, 1))
    end

end


%% Prepare EEG datasets

if strcmp(par.experiment, 'newborn')

    % load channel coordinates file
    load(['/Users/barbero/Documents/Projects/Rhototype_2ioi_newborns/' ...
        'admin/Chan128.mat']);

    % Assign electrodes locations
    for ich = 1:length(chan)

        % Find matching label from lwdata file and channel mat file
        match_label = [];
    
        for i = 1:length(chan)
    
            match_label(i) = strcmp(header.chanlocs(ich).labels, chan(i).labels);
        end
    
        index_ch = find(match_label);
        channel_keep(ich) = index_ch;
    
        % Assignment fields
        header.chanlocs(ich).topo_enabled = 1;
        header.chanlocs(ich).theta = chan(index_ch).theta;
        header.chanlocs(ich).radius = chan(index_ch).radius;
        header.chanlocs(ich).sph_theta = chan(index_ch).sph_theta;
        header.chanlocs(ich).sph_phi = chan(index_ch).sph_phi;
        header.chanlocs(ich).X = chan(index_ch).X;
        header.chanlocs(ich).Y = chan(index_ch).Y;
        header.chanlocs(ich).Z = chan(index_ch).Z;
    
    end

    % keep only EEG channels
    data = data(:,channel_keep,:,:,:,:);
    header.datasize = size(data);
    header.chanlocs = header.chanlocs(channel_keep);


else

    % set channels labels 
    [header, data] = fix_channel_labels(header, data, ...
                                        'chan_labels.csv', ...
                                        par.lw6_path, ...
                                        'keep_eeg_only', false); 

end


%% save 

fpath = fullfile(par.deriv_path, 'merged', sub_num2str(sub)); 
mkdir(fpath); 

header.name = sprintf('sub-%03d_task-EEG', sub);

tic
CLW_save(fpath, header, data); 
fprintf('lw6 save took %.1f seconds. \n', toc); 

