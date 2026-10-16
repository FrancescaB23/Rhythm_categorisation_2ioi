
clear

load_path = '/Users/tomaslenc/projects_backed_up/Rhytotype/data_eeg_2ioi/derivatives/preproc_francesca'; 

subjects = [6 : 23]; 

par = get_par('experiment', '2ioi'); 
 

for i_sub=1:length(subjects)
    
    sub = subjects(i_sub); 
    sub_str = sub_num2str(sub); 
    
    % load and merge time-domain data
    c_ep = 1; 
    events = []; 
    data_all = {}; 
    for i_cond=1:13
        
        % load data 
        fpath = fullfile(par.deriv_path, 'preproc_francesca'); 
        fname =  sprintf('* 99%02d * sub-%03d_task-EEG.lw6', i_cond, sub); 
        d = dir(fullfile(fpath, fname)); 

        [header, data] = CLW_load(fullfile(d.folder, d.name));

        for i_event=1:header.datasize(1)
            events(c_ep).code = num2str(i_cond); 
            events(c_ep).latency = 0; 
            events(c_ep).epoch = c_ep; 
            c_ep = c_ep + 1; 
        end
        
        data_all{i_cond} = data; 
    end
    
    data_ep = cat(1, data_all{:}); 
    header_ep = header; 
    header_ep.datasize = size(data_ep); 
    header_ep.events = events; 
    
    n_cond = length(unique({header_ep.events.code})); 
    assert(n_cond == 13); 

    % save 
    fpath = fullfile(par.deriv_path, 'preproc', sub_str); 
    mkdir(fpath); 
    
    fname = sprintf('sub-%03d_task-EEG_filt_ep_chanint_icfilt', sub); 
    
    header_ep.name = fname; 

    CLW_save(fpath, header_ep, data_ep);
    
end
