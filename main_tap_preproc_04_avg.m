function main_tap_preproc_04_avg(par, sub)
% function to prepare continuous tapping force for further frequency-RSA
% analysis:
% - import data
% - select only orhtogonal tapping force Fz
% - average trials across conditions
% - save data


%%
% load data
fpath = fullfile(par.deriv_path, 'tapping', sub_num2str(sub));
fname = sprintf('sub-%03d_task-TAP_response-force.lw6', sub);
d = dir(fullfile(fpath, fname));

[header_tap, data_tap] = CLW_load(fullfile(d.folder, d.name));


% keep only Fz channel
[header_tap, data_tap] = ...
    RLW_arrange_channels(header_tap, data_tap, {'Fz'});


% average trials per condition
[header_avg, data_avg] = average_trials_per_cond(header_tap, data_tap, par, 'concatenate_epochs', false);

header_avg.name = [header_avg.name, '_avg'];


% save file
CLW_save(fpath, header_avg, data_avg);

end