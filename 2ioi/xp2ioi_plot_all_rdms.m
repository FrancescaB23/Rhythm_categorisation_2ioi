% plot a summary firgure with all matrices for a subset of subjects (either
% based on some criterion of manually defined subjects)

clear 

par = get_par('experiment', '2ioi'); 
subjects = num2cell([6 : 23]); 

n_harmonics = 3; 
harmonics = [1 : n_harmonics];     
harm_str = vec_to_str(harmonics, 'format','%d', 'sep','-'); 


%% eeg lf

response = 'lf'; 
roi_name = 'front';         
avg_chans = true; 

concat_ep = false; 
feat_snr = false; 

if feat_snr 
    snr_str = sprintf('%d-%d', par.snr_bins(1), par.snr_bins(2)); 
else
    snr_str = 'none'; 
end

save_path = fullfile(par.figures_path,...
                         sprintf('response-%s', response), ...
                         sprintf('harm-%s', harm_str), ...
                         sprintf('concatEp-%s_snr-%s', jsonencode(concat_ep), snr_str), ...
                         sprintf('roi-%s_avgChans-%s', roi_name, jsonencode(avg_chans)) ...
                         ); 

fname_save = fullfile(save_path, ...
        sprintf('response-%s_harm-%s_concatEp-%s_snr-%s_roi-%s_avgChans-%s_data.mat', ...
                         response, ...
                         harm_str, ...
                         jsonencode(concat_ep), ...
                         snr_str, ...
                         roi_name, ...
                         jsonencode(avg_chans) ...
                         )); 
tmp = load(fname_save); 
res_lf = tmp.res; 



%%

% idx_to_plot = [res_lf.p_perm] < 0.05; 
idx_to_plot = [1 : length(res_lf)]; 

f = figure('color', 'white', 'position', [620 2 160 1003]); 
pnl = panel(f); 

pnl.pack('h', 1); 

plot_all_rdms(pnl(1), ...
    [res_lf(idx_to_plot).sub], ...
    {res_lf(idx_to_plot).rdm}, ...
    {res_lf(idx_to_plot).best_model}, ...
    [res_lf(idx_to_plot).p_perm], ...
    [res_lf(idx_to_plot).z_snr] ...
    )


pnl.de.margin = [0, 0, 0, 8]; 
pnl(1).marginright = 10; 
pnl.margin = [10, 2, 10, 15]; 


%%
% 
% fname = fullfile(par.figures_path, 'hf-significant_subjects.png'); 
% 
% fname = fullfile(par.figures_path, 'all_subjects.png'); 
% 
% print(f, '-dpng', '-painters', '-r300', fname); 
% 
