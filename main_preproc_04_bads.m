% This script helps to do data cleaning, to visually identify bad channels
% and trials that can be subsequently interpolated or removed, respectively.

clear 
clc


par = get_par('experiment', 'newborn'); % 'newborn' '2ioi_conga' 

sub = 28
% 'newborn' - subjects = [1:3 5:14 16:25 27:29 31:42];
% '2ioi_conga' - subjects = [2:20]; % [2:7 9:15 17:20]


%% manual check for bad chans and trials 

wd = pwd;

% use letswave 7 to open continuous data viewer and check for bad
% trials/channels -> manually write them down 
import_lw(7); 
cd(fullfile(par.deriv_path, 'preproc', sub_num2str(sub))); 
letswave 

warning('set X-range in the viewer to 0 - %.1fs', par.trial_dur); 


%%  
% switch back to letswave 6 when done 
cd(wd); 
import_lw(6);

        
%% generate template table for bad trials and channels 

% This chunk of code will look whether the text file with info about bad
% chanenls and trials exists for the current subject. If not, it will
% generate a new empty one in the proper location and with a proper name.
% Once this is done, open the files manually and enter the bad channels and
% trials you observed during visual inspection in the previous step.

% prepare info about bad channels
fname_bad_chans = fullfile(par.deriv_path, 'preproc', sub_num2str(sub),...
                 sprintf('%s_badChans.csv', sub_num2str(sub))); 
             
if ~exist(fname_bad_chans)
    tbl_bad_chans = cell2table(cell(0, 1), 'VariableNames', {'channel'}); 
    writetable(tbl_bad_chans, fname_bad_chans); 
end

% prepare info about bad trials
fname_bad_trials = fullfile(par.deriv_path, 'preproc', sub_num2str(sub), ...
                 sprintf('%s_badTrials.csv', sub_num2str(sub))); 
             
if ~exist(fname_bad_trials)
    tbl_bad_trials = cell2table(cell(0, 1), 'VariableNames', {'trial'}); 
    writetable(tbl_bad_trials, fname_bad_trials); 
end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% EDIT THE TEXT FILE MANUALLY TO ADD/REMOVE BADS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%%%
% reject/fix bads

% load the info about bad channels and trials
tbl_bad_chans = readtable(fname_bad_chans); 
tbl_bad_trials = readtable(fname_bad_trials); 

bad_chans = tbl_bad_chans.channel;
bad_epochs = tbl_bad_trials.trial; 

% print to console for sanity check 
fprintf('\n\n---------------------------------------------------------\n'); 
bad_chans
bad_epochs
fprintf('---------------------------------------------------------\n'); 


%%
% load EEG data

if strcmp(par.experiment, 'newborn')

    fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
    fname =  sprintf('sub-%03d_task-EEG_filt_ep_AB-100', sub); % AB cleaned files, threshold 100 microV
    [header_ep, data_ep] = CLW_load(fullfile(fpath, fname));

else

    fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
    fname =  sprintf('sub-%03d_task-EEG_filt_ep', sub); 
    [header_ep, data_ep] = CLW_load(fullfile(fpath, fname));

end


%% Make sure there are enough repetitions per condition (at least one)

% if it is not the case, reconsider rejected trials

header_check_trials = header_ep;
idx_keep_trials = setdiff([1:header_ep.datasize(1)], bad_epochs);
header_check_trials.events = header_check_trials.events(idx_keep_trials);
cond_code = unique({header_ep.events.code});

for i = 1:length(cond_code)

    rep_cond(i) = sum(strcmp({header_check_trials.events.code}, ...
                                                    cond_code(i))); 
end

disp_trial_info = table(string(cond_code).', rep_cond', ...
    'VariableNames', ["condition","repetition"]);
disp(disp_trial_info);

cond_keep_code = unique({header_check_trials.events.code});

% If not enough trial for at least one condition, show message
if length(cond_keep_code) < 10

    disp(['Not enough trials per condition, reject less! ' ...
        'Missing trials for conditions:']);

    cond_miss = setdiff(cond_code, cond_keep_code);
    disp(cond_miss);

    % show trial index for missing conditions, so that at least one trial
    % can be reintegrated in the data
    for j = 1:length(cond_miss)

        disp('epoch index per condition:')
        disp(cond_miss(j))
        disp(find(ismember({header_ep.events.code}, cond_miss(j))))
    end
  
end

% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
% EDIT THE TEXT FILE MANUALLY FOR BAD TRIALS AND RE-RUN LOAD OF BAD TRIALS
% !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


%%
% reject bad trials 
[header_ep, data_ep] = RLW_arrange_epochs(header_ep, data_ep, ...
                            setdiff([1:header_ep.datasize(1)], bad_epochs)); 

% interpolate bad channels 
if ~isempty(bad_chans)
    [header_ep, data_ep] = interpolate_bad_chans(header_ep, data_ep, ...
        bad_chans, par.n_chans_interp); 
    
    bad_chans_str = join(bad_chans, '-');

end

% save
header_ep.name = sprintf('%s_chanint', header_ep.name);
CLW_save(fpath, header_ep, data_ep);


size(data_ep)