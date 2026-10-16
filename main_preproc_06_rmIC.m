% This script can be used to visually identify bad IC components and remove
% them from the data

clear 
clc

par = get_par('experiment', '2ioi_conga'); 

sub = 20


%% visual inspecion of ICs 

% open letswave6, and manually inspect the ICs in the
% "sub-$_task-EEG_filt_ep_chanint" dataset, and write down which ones are
% bad

wd = pwd; 
fpath = fullfile(par.deriv_path, 'preproc', sub_num2str(sub)); 
cd(fpath); 
letswave 


%% change dir back to code

cd(wd); 


%% save/load bad ICs 
% Same as for bad channels and tials, this chunk of code will look whether
% the text file with data on bad ICs exists. If not it will create an empty
% one that you can edit manually.

fname = fullfile(par.deriv_path, 'preproc', sub_num2str(sub), ...
                 sprintf('%s_badICs.csv', sub_num2str(sub))); 
if ~exist(fname)
    tbl_bad_ICs = cell2table(cell(0, 1), 'VariableNames', {'IC'}); 
    writetable(tbl_bad_ICs, fname); 
end

% EDIT THE TEXT FILE MANUALLY

%

% load the bad IC info from file 
tbl_bad_ICs = readtable(fname); 
bad_ICs = tbl_bad_ICs.IC;

% print to console for sanity check 
fprintf('\n\n---------------------------------------------------------\n'); 
bad_ICs
fprintf('---------------------------------------------------------\n'); 


%%
% load the ICA matrix 
clear matrix_ICA
fname = sprintf('sub-%03d_task-EEG_icaMatrix', sub); 
load(fullfile(fpath, fname));

% load data 
fname =  sprintf('sub-%03d_task-EEG_filt_ep_chanint', sub); 
[header_ep, data_ep] = CLW_load(fullfile(fpath, fname));

% remove bad ICs from data 
n_ic = size(matrix_ICA.ica_mm, 2); 
ic2keep = setdiff([1 : n_ic], bad_ICs); 

[header_ep, data_ep] = RLW_ICA_apply_filter(header_ep, data_ep, ...
                                      matrix_ICA.ica_mm, ...
                                      matrix_ICA.ica_um, ...
                                      ic2keep);     

% save clean data 
header_ep.name = [header_ep.name, sprintf('_icfilt(%s)', ...
                                          vec_to_str(bad_ICs, 'sep',','))];

CLW_save(fpath, header_ep, data_ep); 
