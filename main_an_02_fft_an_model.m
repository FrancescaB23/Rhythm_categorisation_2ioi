% Run FFT on the auditory-nerve model output. 

clear

par = get_par('experiment', '2ioi_conga'); 

%% 

% prepare output directories
save_path = fullfile(par.deriv_path, 'fft', 'response-urear-an'); 
mkdir(save_path); 


%% load auditory nerve responses in a letswave format

load_path = fullfile(par.deriv_path, 'cochlear_model'); 

an = cell(1, par.n_cond); 
events = []; 

for i_cond=1:par.n_cond
    
    fname = sprintf('response-AN_cond-%02d.mat', i_cond); 
    
    res = load(fullfile(load_path, fname)); 
    an{i_cond} = sum(res.AN.an_sout, 1); 
    
    events(i_cond).code = num2str(i_cond); 
    events(i_cond).epoch = i_cond; 
    events(i_cond).latency = 0; 
    
end

data(:, 1, 1, 1, 1, :) = cat(1, an{:}); 

header = make_lw_header('data', data, 'fs', res.AN.Fs_an_sout, 'events', events); 
    

%% compute FFT and save

fname = sprintf('response-urear-an'); 

fmax = (1/header.xstep) / 2;

% magnitude FFT
[header_mX, data_mX] = RLW_FFT(header, data);

[header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax);

data_mX(:, :, :, :, :, 1) = 0;

header_mX.name = 'response-urear-an_snr-none_mX';

CLW_save(save_path, header_mX, data_mX); 
                           

% complex FFT
[header_X, data_X] = RLW_FFT(header, data, 'output', 'complex');

[header_X, data_X] = crop_lw_data(header_X, data_X, fmax);

data_X(:, :, :, :, :, 1) = 0;

header_X.name = 'response-urear-an_snr-none_X';

CLW_save(save_path, header_X, data_X); 