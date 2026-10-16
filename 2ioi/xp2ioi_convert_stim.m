
clear 

subjects = 1;%[6 : 23]; 

par = get_par('experiment', '2ioi_vision'); %'2ioi'


%%

% load stimulus info 
% load('/Users/tomaslenc/projects_backed_up/Rhytotype/data_eeg_2ioi/stimuli/medium-tempo_67.00-33.00_tone300_050ms_22.50s_0.85int/medium-tempo_67.00-33.00_tone300_050ms_22.50s_0.85int.mat'); 
% 2ioi_vision
load(['/Users/barbero/Documents/GitKraken/Projects/' ...
    'Rhototype_2ioi_vision/stimulation/eeg/sequence_for_visual/' ...
    'sequences_drop_volume/drop_0.85int_13cond_240ms_22.50s/' ...
    'medium-tempo_67.00-33.00_tone300_240ms_22.50s_0.85int/' ...
    'medium-tempo_67.00-33.00_tone300_240ms_22.50s_0.85int.mat']);

ratios = cfg.perc / 100; 
ratio_labs = get_ratio_labs(ratios(:, 1)); 

stim_all = num2cell(cfg.drop.sequence_zero_drop, 2); 

stim = []; 
stim.s = stim_all; 
stim.fs = cfg.Fs; 
stim.ratios = ratios; 
stim.ratio_labs = ratio_labs; 

stim.cycle_dur = cfg.patternPeriod; 
stim.trial_dur = cfg.seqDuration; 
stim.n_cond = length(ratios); 

%%

for i_sub=1:length(subjects)
    
    sub = subjects(i_sub); 
    sub_str = sub_num2str(sub); 

    fpath = fullfile(par.deriv_path, 'stimuli', sub_str); 
    fname = sprintf('sub-%03d_stim', sub);

    mkdir(fpath); 
    save(fullfile(fpath, fname), 'stim')
    
end