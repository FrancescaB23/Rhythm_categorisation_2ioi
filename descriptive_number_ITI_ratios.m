% script to retrieve the number of valid tap ratios 

clc
clear
close all


%% parameters
par = get_par('experiment', '2ioi_conga');
subjects = num2cell([2:7 9:15 17:20]);


%%

for i_sub = 1:length(subjects)

    sub = subjects{i_sub};
    sub_str = sub_num2str(sub);

    load_path = fullfile(par.deriv_path, 'tap_onsets', sub_str); 

    fname = sprintf('sub-%03d_response-tap-ratios.mat', sub); 

    res = load(fullfile(load_path, fname));

    n_valid_tap_ratios(i_sub, :) = res.n_valid_tap_ratios;

end

%% 

fprintf('\nmean number valid ITI ratios per condition: %.4f', ...
    mean(mean(n_valid_tap_ratios)))
fprintf('\nminimum number valid ITI ratios per condition: %.4f', ...
    max(max(n_valid_tap_ratios)))
fprintf('\nmaximum number valid ITI ratios per condition: %.4f\n', ...
    min(min(n_valid_tap_ratios)))

