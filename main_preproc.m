% Script to call some preprocessing functions to loop across participants

clc
clear
close all


%% parameters
par = get_par('experiment', '2ioi_conga'); % 'newborn' '2ioi_conga'

% specify participants depending on experiment
switch par.experiment

    case 'newborn'
        subjects = [1:3 5:14 16:25 27:29 31:42];

    case '2ioi_conga'
        subjects = [2:7 9:15 17:20];
end


%%
for i_sub = 1:length(subjects)

    % eeg

%     % filter
%     disp(sprintf('filter - group-%s sub-%03d', ...
%         par.experiment, subjects(i_sub)));
%     main_preproc_02_filter(par, subjects(i_sub));
% 
% 
%     % apply Artefact Blocking algorithm on newborn data
%     if strcmp(par.experiment, 'newborn')
%         disp(sprintf('AB algorithm - group-%s sub-%03d', ...
%             par.experiment, subjects(i_sub)));
%         main_preproc_03_AB(par, subjects(i_sub));
%     end


%     % perform ICA
%     disp(sprintf('ICA - group-%s sub-%03d', par.experiment, subjects(i_sub)));
%     main_preproc_05_ica(par, subjects(i_sub));


%     % (remove more external electrodes - newborns only), segment,
%     % downsample, rereference to mastoids, and average trials per condition
%     disp(sprintf(['rereference, average trials per condition, and more - ' ...
%         'group-%s sub-%03d'], par.experiment, subjects(i_sub)));
%     main_preproc_07_rr_avg(par, subjects(i_sub));

%     % erp adults (all chunks included)
%     disp(sprintf('EEG erp - group-%s sub-%03d', par.experiment, subjects(i_sub)));
%     main_preproc_11_erp(par, subjects(i_sub));


    % tapping

%     % average continuous tapping force signal
%     disp(sprintf(['tapping: average trials per condition - ' ...
%         'group-%s sub-%03d'], par.experiment, subjects(i_sub)));
%     main_tap_preproc_04_avg(par, subjects(i_sub));

    % erp (all chunks included
    disp(sprintf('TAP erp - group-%s sub-%03d', par.experiment, subjects(i_sub)));
    main_tap_preproc_05_erp(par, subjects(i_sub));

end