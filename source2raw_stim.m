% This script does kind of source2raw conversion for the stimulus audio
% data. We will assume that the stimuls audio was identical for all
% subjects in the experiment.

clear 

par = get_par('experiment', 'newborn'); %  newborn 2ioi_conga 


% the most convenient way is to collect all the parameters and audio
% stimulus signals into a single .mat file. The best file
% to take for this is the .mat stimulus file for the last condition...
% let's find it and load it

if strcmp(par.experiment,'newborn')

    load(['/Users/barbero/Documents/Projects/Rhototype_2ioi_newborns/' ...
        'newborns/pilot_testing/stimuli/dropInt-0.85_nCond-10_' ...
        'sample-conga_eventDur-200ms_patDur-750ms_seqDur-30000ms/' ...
        'iCond-10_ratios-0.667-0.333/' ...
        'iCond-10_ratios-0.667-0.333.mat']);

elseif strcmp(par.experiment, '2ioi_conga')

    % pick one subject 
    sub = 2; 
    
    d_folders = dir(fullfile(par.raw_path, sub_num2str(sub), 'stimuli','iCond*')); 
    
    cond_nums = cellfun(@(x) str2num(x{:}), ...
        regexp({d_folders.name}, '(?<=iCond-)\d+', 'match')); 
    
    n_cond = max(cond_nums); 
    
    d = dir(fullfile(par.raw_path, sub_num2str(sub), 'stimuli',...
        sprintf('iCond-%d_*', n_cond))); 
    
    assert(length(d) == 1, ...
        'could not find last condition stimulus in raw/stimuli'); 
    
    d = dir(fullfile(d.folder, d.name, '*.mat')); 
         
    load(fullfile(d.folder, d.name)); 

end

% prepare new structure that will have all the relevant information
% available (audio data, parameters, etc.)
if ~isfield(cfg, 'ratios')
    cfg.ratios = cfg.perc / 100;
end
ratios = cfg.ratios; 
ratio_labs = get_ratio_labs(ratios(:, 1)); 
if isfield(cfg.drop, 'sequence_zero_drop')
    stim_all = num2cell(cfg.drop.sequence_zero_drop, 2); 
else
    stim_all = 0;
end

stim = []; 
stim.s = stim_all; 
if isfield(cfg, 'fs')
    stim.fs = cfg.fs; 
else 
    stim.fs = cfg.Fs;
end
stim.ratios = ratios; 
stim.ratio_labs = ratio_labs; 


if isfield(cfg, 'pattern_period')
    stim.cycle_dur = cfg.pattern_period; 
else
    stim.cycle_dur = cfg.patternPeriod; 
end
if isfield(cfg, 'seq_duration')
    stim.trial_dur = cfg.seq_duration; 
else
    stim.trial_dur = cfg.seqDuration; 
end
if isfield(cfg, 'n_cond')
    stim.n_cond = cfg.n_cond; 
else
    stim.n_cond = cfg.nCond;
end

% save ti to derivatives
fpath = fullfile(par.deriv_path, 'stimuli'); 
mkdir(fpath); 
save(fullfile(fpath, 'stim'), 'stim')