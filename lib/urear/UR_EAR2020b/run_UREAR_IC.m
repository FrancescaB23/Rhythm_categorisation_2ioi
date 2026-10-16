function [IC, cfg] = run_UREAR_IC(cfg, AN, BMFs)
            
cfg.BMFs = BMFs; 

% allocate IC
% (there is 6-7 ms added in the output so we can allocate more...cut it
% within the loop before assigning)
numCF                   = length(cfg.CFs); 
numBMF                  = length(cfg.BMFs);    
N                       = round(size(AN.an_sout, 2) / AN.Fs_an_sout * cfg.Fs_decim);  
BE_sout_population      = nan(numBMF, numCF, N); 
BS_sout_population      = nan(numBMF, numCF, N); 

for bmfi=1:numBMF

    BMF = cfg.BMFs(bmfi); 
    fprintf('processing IC (iter %d) with BMF = %g Hz \n',bmfi,BMF); 

    for cfi=1:numCF

        switch cfg.Which_IC

            case 1 % Monaural SFIE

                [ic_sout_BE, ic_sout_BS, cn_sout_contra] = SFIE_BE_BS_BMF(...
                    AN.an_sout_spikes(cfi,:), BMF, AN.Fs_an_sout_spikes);
                
                BE_sout_population_tmp = resample(ic_sout_BE, cfg.Fs_decim, AN.Fs_an_sout_spikes);
                BE_sout_population(bmfi,cfi,:) = BE_sout_population_tmp(1:N);

                BS_sout_population_tmp = resample(ic_sout_BS, cfg.Fs_decim, AN.Fs_an_sout_spikes);
                BS_sout_population(bmfi,cfi,:) = BS_sout_population_tmp(1:N);

            case 2 % Monaural Simple Filter

                % Now, call NEW unitgain BP filter to simulate bandpass IC cell with all BMFs.
                ic_sout_BE = unitgain_bpFilter(AN.an_sout_spikes(cfi,:), BMF, AN.Fs_an_sout_spikes); 

                BE_sout_population_tmp = resample(ic_sout_BE, cfg.Fs_decim, AN.Fs_an_sout_spikes);
                BE_sout_population(bmfi,cfi,:) = BE_sout_population_tmp(1:N);
        end

    end

end

IC           = []; 
IC.BE_sout   = BE_sout_population;
IC.Fs        = cfg.Fs_decim; 



