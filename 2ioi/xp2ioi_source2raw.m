
clear

par = get_par('experiment', '2ioi'); 

subjects = [6 : 23]; 

source_path = fullfile(par.data_path, 'source'); 


for i_sub=1:length(subjects)
    
    sub = subjects(i_sub); 
    sub_str = sub_num2str(sub); 
    
    % EEG (ignore TAP, Fra did it)
    d = dir(sprintf('%s/eeg/%s/*EEG*.bdf', source_path, sub_str)); 
    
    save_path = fullfile(par.raw_path, sub_str, 'eeg'); 
    mkdir(save_path); 

    for i_file=1:length(d)
        copyfile(fullfile(d(i_file).folder, d(i_file).name), ...
                 fullfile(save_path, d(i_file).name)); 
    end
    
    % log
    d = dir(sprintf('%s/logfile/%s/*EEG*.tsv', source_path, sub_str)); 
    
    save_path = fullfile(par.raw_path, sub_str, 'log'); 
    mkdir(save_path); 

    for i_file=1:length(d)
        copyfile(fullfile(d(i_file).folder, d(i_file).name), ...
                 fullfile(save_path, d(i_file).name)); 
    end

end

copyfile(fullfile(source_path, 'README.rtf'), ...
         fullfile(par.raw_path, 'README.rtf')); 