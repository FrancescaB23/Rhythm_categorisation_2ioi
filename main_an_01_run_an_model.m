% Run auditory nerve model

clear

par = get_par('experiment', '2ioi_conga'); 

%% model parameters

cfg = []; 

% sampling frequencies
% Model sampling rate (must be 100k, 200k or 500k for AN model):
cfg.Fs_wav_input = 44100; % oroginal input sampling rate
cfg.Fs           = 100e3; % samples/sec
cfg.Fs_decim     = 2000;  % resample rate for plots

cfg.spl = 73.2;  % dB SPL of the stimulus

% AN fiber type. (1 = low SR, 2 = medium SR, 3 = high SR) (ONLY for Zilany 2014 model!!!)
cfg.fiberType = 3; 

% Set cfg.implnt = 0 => approximate model, 1 => exact power law
% implementation (See Zilany et al., 2009)
cfg.implnt = 0;

% number of repetitions
cfg.nrep = 1;

% Set cfg.noiseType to 0 for fixed fGn or 1 for variable fGn - this is the
% 'noise' associated with spontaneous activity of AN fibers - see Zilany
% et al., 2009. 0 lets you "freeze" it.
cfg.noiseType = 1;

% 1=cat; 2=human AN model parameters (with Shera tuning sharpness)
cfg.species = 2; 

% Number of [LSR, MSR, HSR] fibers at each CF in a healthy AN (only for Bruce 2018 model)
cfg.n_fibres = 51; % 51 (total number of fibres, deprecated)
cfg.numsponts = round([cfg.n_fibres*0.16, cfg.n_fibres*0.23, cfg.n_fibres*0.61]); % (based on Liberman 1978)

% set range and resolution ofscfg.CFs here
cfg.CF_num = 128; % 128
cfg.CF_range = [130, 3000]; 
cfg.CFs = logspace(log10(cfg.CF_range(1)), log10(cfg.CF_range(2)), cfg.CF_num); 
cfg.CFs([1 end]) = cfg.CF_range; % force end points to be exact.

% PSTH parameters (only for Bruce 2018)
cfg.psthbinwidth_mr = 100e-6; % mean-rate binwidth in seconds;
cfg.windur_ft       = 32;
cfg.smw_ft          = hamming(cfg.windur_ft);
cfg.windur_mr       = 128;
cfg.smw_mr          = hamming(cfg.windur_mr);
cfg.psthbins        = round(cfg.psthbinwidth_mr * cfg.Fs);  % number of psth_ft bins per psth bin

% audiogram (normal hearing)
cfg.ag_fs = [125, 250, 500, 1000, 2000, 4000, 8000]; 
cfg.ag_dbloss = [0 0 0 0 0 0 0];
cfg.dbloss = interp1(cfg.ag_fs, cfg.ag_dbloss, cfg.CFs, 'linear', 'extrap');
[cfg.cohc_vals, cfg.cihc_vals] = fitaudiogram2(cfg.CFs, cfg.dbloss, cfg.species);

% For a very low CF, a 0 may be returned by Bruce et al. fit
% audiogram, but this is a bad default.  Set it to 1 here.
if cfg.cohc_vals(1) == 0
    cfg.cohc_vals(1) = 1;
end
if cfg.cihc_vals(1) == 0
    cfg.cihc_vals(1) = 1;
end

% define which AN model will be used: 
%     1: Zilany 2014
%     2: Bruce 2018
cfg.Which_AN = 2; 


%%

load(fullfile(par.deriv_path, 'stimuli', 'stim.mat')); 

save_path = fullfile(par.deriv_path, 'cochlear_model'); 

if ~isdir(save_path)
    mkdir(save_path); 
end

for i_cond=1:length(stim.s)

    s = stim.s{i_cond}; 
    fs = stim.fs; 

    assert(fs == cfg.Fs_wav_input);
   
    % run
    [AN, par_model] = run_UREAR_AN(ensure_row(s), cfg); 

    % save
    fname = sprintf('response-AN_cond-%02d.mat', i_cond); 

    save(fullfile(save_path, fname), 's', 'fs', 'AN', 'par_model', '-v7.3')

end