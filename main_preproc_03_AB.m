function main_preproc_03_AB(par, sub)
% function to apply Artefact Blocking algorithm to newborn data


%% 
% load
fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 

fname = sprintf('sub-%03d_task-EEG_filt_ep', sub); 

[header, data] = CLW_load(fullfile(fpath, fname)); 


% Apply artefact blocking (AB) algorithm 

addpath(genpath(par.ab_path)); 

data_tmp = squeeze(permute(data,[6,5,4,3,1,2])); 
data_concat = reshape(data_tmp, [], size(data,2))'; 

cfg = [];
cfg.Approach = 'Window'; 
cfg.Threshold = 100; % the threshold for detecting artifact
cfg.Fs = 1/header.xstep; % sampling frequency
cfg.WindowSize = 10; % unit in second
cfg.InData = data_concat;  % the EEG data

cfg = Run_AB(cfg);

n_epochs = size(data,1); 
ab_data = reshape(cfg.OutData', header.datasize(end), n_epochs, size(data,2));
ab_data = permute(ab_data, [2,3,1]); 
out_data = []; 
out_data(:,:,1,1,1,:) = ab_data;

out_header = header;
out_header.name = sprintf('sub-%03d_task-EEG_filt_ep_AB-%d', sub, cfg.Threshold); 

CLW_save(fpath, out_header, out_data); 


end