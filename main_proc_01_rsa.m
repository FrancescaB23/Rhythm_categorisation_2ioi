% script to perform RSA from time domain signal

% NOTE: correlation with only 2 values can result in either +1 or -1. So
% beware that RDM obtained with magnitudes of first two harmonics will be
% full of +1 and -1 (sometimes only +1 and then you'll see a 2"flat"
% matrix).

clc
close all
clear


%% parameters

par = get_par('experiment', '2ioi_conga'); % 'newborn' '2ioi_conga'

feat_names = {'mX'}; % {'mX', 'X'}; %

par_rsa = [];

par_rsa.harmonics = [1:6];

% feat_name_stim is specified in the script

flip_rdm = false;

rerun_stats = false;
rerun_figs = false;


%%
switch par.experiment

    case 'newborn'

        subjects = {'grand-19'}; %[{'grand-19'}, {'grand-38-zSNR'}];
        %         subjects = num2cell([1 5 7 8 9 10 12 17 18 22 23 27 28 31 33 35 36 37 38]);

    case '2ioi_conga'

        subjects = [num2cell([2:7 9:15 17:20]), {'grand'}];
end


%%

for i_feat=1:length(feat_names)

    feat_name_resp = feat_names{i_feat};

    feat_name_stim = 'ratio';
%         feat_name_stim = feat_names{i_feat};
%         feat_name_stim = 'cochlear';


    for i_sub=1:length(subjects)

        sub = subjects{i_sub};


        %% EEG FFT
%         for feat_snr = false
% 
%             par_rsa.feat_snr = feat_snr;
% 
%             % --> frontal ROI
%             par_rsa.avg_chans = true;
%             par_rsa.roi_name = 'front';
% 
%             rsa_lf_time(sub, par_rsa, ...
%                 'rerun_stats', rerun_stats, ...
%                 'rerun_figs', rerun_figs, ...
%                 'flip_rdm', flip_rdm, ...
%                 'feat_resp', feat_name_resp, ...
%                 'feat_stim', feat_name_stim, ...
%                 'par', par);
%         end


        %% EEG ERP
%         % --> frontal ROI
%         par_rsa.avg_chans = true;
%         par_rsa.roi_name = 'front';
% 
%         rsa_erp(sub, par_rsa, ...
%             'rerun_stats', rerun_stats, ...
%             'rerun_figs', rerun_figs, ...
%             'flip_rdm', flip_rdm, ...
%             'feat_stim', feat_name_stim, ...
%             'par', par);


        %% TAP RATIO
%         rsa_tap_ratios(sub, ...
%             'rerun_stats', rerun_stats, ...
%             'rerun_figs', rerun_figs, ...
%             'flip_rdm', flip_rdm, ...
%             'feat_stim', feat_name_stim, ...
%             'par', par);


        %% TAP FORCE
%         par_rsa.feat_snr = false;
% 
%         rsa_tap_force(sub, par_rsa, ...
%             'rerun_stats', rerun_stats, ...
%             'rerun_figs', rerun_figs, ...
%             'flip_rdm', flip_rdm, ...
%             'feat_resp', feat_name_resp, ...
%             'feat_stim', feat_name_stim, ...
%             'par', par);

            
        %% TAP ERP
        rsa_tap_force_erp(sub, par_rsa, ...
            'rerun_stats', rerun_stats, ... 
            'rerun_figs', rerun_figs, ...
            'flip_rdm', flip_rdm, ... 
            'feat_stim', feat_name_stim, ...
            'par', par);
        %%%%
    end
end