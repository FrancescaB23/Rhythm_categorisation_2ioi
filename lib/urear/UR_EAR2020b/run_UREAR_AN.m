function [AN,cfg] = run_UREAR_AN(s, cfg)

Pref = 20e-6; % reference pressure in pascals

% resample to sampling rate required for AN model, and store in
s = resample(s, cfg.Fs, cfg.Fs_wav_input);

% scale to desired RMS level
s = s * (Pref * 10.^(cfg.spl/20) / rms(s));

% stimulus length
N               = length(s);
stimDur         = N/cfg.Fs;
stimDurPadded   = stimDur + 0.04;
T               = 1/cfg.Fs;

% check for NaN in the stimulus (otherwise Matlab will pass these to mex
% filex and crash hard)
if any(isnan(s))
    error('NaN in the stimulus')
end

%% model AN stage

% Run through this once, because this will be the same for any
% IC unit. 
% After the AN is done, go to another, separate loop across IC units. 

% initialize variables
numCF = length(cfg.CFs);
cfis = 1:numCF; 

% allocate AN
if cfg.Which_AN==1

    N                   = round(stimDurPadded * cfg.Fs_decim);  
    an_sout_plot        = nan(numCF, N); % this will be for plotting
    an_sout             = zeros(numCF, stimDurPadded*cfg.Fs); % this will be input to IC unit

elseif cfg.Which_AN==2

    [sponts,tabss,trels] = generateANpopulation(cfg.CF_num, cfg.numsponts);

    CFlp=1; spontlp=1; CF = cfg.CFs(CFlp);

    sponts_concat       = [sponts.LS(CFlp,1:cfg.numsponts(1)) sponts.MS(CFlp,1:cfg.numsponts(2)) sponts.HS(CFlp,1:cfg.numsponts(3))];
    tabss_concat        = [tabss.LS(CFlp,1:cfg.numsponts(1))  tabss.MS(CFlp,1:cfg.numsponts(2))  tabss.HS(CFlp,1:cfg.numsponts(3))];
    trels_concat        = [trels.LS(CFlp,1:cfg.numsponts(1))  trels.MS(CFlp,1:cfg.numsponts(2))  trels.HS(CFlp,1:cfg.numsponts(3))];

    cohc = cfg.cohc_vals(CFlp);
    cihc = cfg.cihc_vals(CFlp);

    [vihc]              = model_IHC_BEZ2018(s,CF,cfg.nrep,T,stimDurPadded,cohc,cihc,cfg.species); 

    spont               = sponts_concat(spontlp);
    tabs                = tabss_concat(spontlp);
    trel                = trels_concat(spontlp);    

    [psth_ft]           = model_Synapse_BEZ2018(vihc, CF, cfg.nrep, T, cfg.noiseType, cfg.implnt, spont, tabs, trel);

    psth_mr             = sum( reshape( psth_ft, cfg.psthbins, length(psth_ft)/cfg.psthbins ) );

    neurogram_ft        = zeros(numCF, length(psth_ft)); % fuck this must be initialized to zeros, do'nt do NaN cause you're summing!!!
    neurogram_mr        = zeros(numCF, length(psth_mr)); 

    cfg.FsFt            = 1 / 1.6000e-04; 
    cfg.FsMr            = 1 / 0.0064; 

    N                   = ceil(stimDurPadded * cfg.FsMr);  
    an_sout_plot        = zeros(numCF, N); % this will be the mean rate (to save space) 
    an_sout             = zeros(numCF, length(psth_ft)); % this will be input to IC unit

end

% open parpool if not open yet
% if isempty(gcp('nocreate'))
%     parpool(4); 
% end

% do not open parpool and change from parfor to for as it is not working
% locally

% CF loop
% parfor cfi=cfis
for cfi=cfis

    % Get one element of each array.
    CF = cfg.CFs(cfi); % CF in Hz;
    fprintf('processing channel %d with CF = %g Hz \n',cfi,CF); 

    cohc = cfg.cohc_vals(cfi);
    cihc = cfg.cihc_vals(cfi);

    switch cfg.Which_AN

        case 1

            % Using ANModel_2014 (2-step process)
            vihc = model_IHC(s,CF,cfg.nrep,T,stimDurPadded,cohc,cihc,cfg.species);

            % an_sout is the auditory-nerve synapse output - a rate vs. time
            % function that could be used to drive a spike generator.
            [an_sout(cfi,:),~,~] = model_Synapse(vihc,CF,cfg.nrep,T,cfg.fiberType,cfg.noiseType,cfg.implnt);

            % Save synapse output waveform into a matrix.
            an_sout_plot(cfi,:) = resample(an_sout(cfi,:), cfg.Fs_decim, cfg.Fs);


        case 2

            % 2018 model
            vihc = model_IHC_BEZ2018(s,CF,cfg.nrep,T,stimDurPadded,cohc,cihc,cfg.species);      


            for spontlp = 1:sum(cfg.numsponts)

                spont = sponts_concat(spontlp);
                tabs = tabss_concat(spontlp);
                trel = trels_concat(spontlp);

                [psth_ft,~,~,~] = model_Synapse_BEZ2018(vihc, CF, cfg.nrep, T, cfg.noiseType, cfg.implnt, spont, tabs, trel);

                if spontlp==1
                    psth = psth_ft; 
                else
                    psth = psth + psth_ft; 
                end

                psth_mr  = sum( reshape( psth_ft, cfg.psthbins, length(psth_ft)/cfg.psthbins ) );

                neurogram_ft(cfi,:) = neurogram_ft(cfi,:) + filter(cfg.smw_ft,1,psth_ft);
                neurogram_mr(cfi,:) = neurogram_mr(cfi,:) + filter(cfg.smw_mr,1,psth_mr);

            end 

            an_sout(cfi,:) = (100000*psth)/sum(cfg.numsponts);

    end

end % end of CF loop


% if 2018 AN model
if cfg.Which_AN==2

    neurogram_ft = neurogram_ft(:, 1:cfg.windur_ft/2:end); % 50% overlap in Hamming window
    t_ft = [...
        0 : ...
        cfg.windur_ft / 2 / cfg.Fs : ...
        (size(neurogram_ft,2)-1) * cfg.windur_ft / 2 / cfg.Fs...
        ]; % time vector for the fine-timing neurogram

    neurogram_mr = neurogram_mr(:, 1:cfg.windur_mr/2:end); % 50% overlap in Hamming window
    t_mr = [...
        0 : ...
        cfg.windur_mr / 2 * cfg.psthbinwidth_mr : ...
        (size(neurogram_mr,2)-1) * cfg.windur_mr / 2 * cfg.psthbinwidth_mr...
        ]; % time vector for the mean-rate neurogram

end


% we need to cut the silence at the end (we put 40 ms)
AN = []; 

if cfg.Which_AN == 1

    maxIdx              = round(stimDur * cfg.Fs_decim); 
    data                = an_sout_plot(:,1:maxIdx);

    AN.an_sout_spikes    = an_sout; 
    AN.an_sout           = data; 
    AN.Fs_an_sout_spikes = cfg.Fs;     
    AN.Fs_an_sout        = cfg.Fs_decim; 

elseif cfg.Which_AN == 2

    maxIdx              = round(stimDur * cfg.FsFt); 
    data                = neurogram_ft(:,1:maxIdx);

    AN.an_sout_spikes    = an_sout; 
    AN.an_sout           = data; 
    AN.Fs_an_sout_spikes = cfg.Fs;     
    AN.Fs_an_sout        = cfg.FsFt;     

end

