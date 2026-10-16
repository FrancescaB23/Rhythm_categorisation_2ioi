function out = get_and_save_ffts(header, data, snr_bins, varargin)

parser = inputParser; 

addParameter(parser, 'fmax', (1/header.xstep) / 2)
addParameter(parser, 'save', true)
addParameter(parser, 'fname', [])
addParameter(parser, 'save_path', [])

parse(parser, varargin{:}); 

fmax = parser.Results.fmax; 
do_save = parser.Results.save; 
fname = parser.Results.fname; 
save_path = parser.Results.save_path; 


if do_save && (isempty(fname) || isempty(save_path))
    error('cannot save witohut save_path and fname provided!');     
end

% ---- magnitude FFT -----

[header_mX, data_mX] = RLW_FFT(header, data); 

[header_mX, data_mX] = crop_lw_data(header_mX, data_mX, fmax); 

data_mX(:, :, :, :, :, 1) = 0; 

header_mX.name = sprintf('%s_snr-none_mX', fname); 

if do_save
    CLW_save(save_path, header_mX, data_mX); 
end

% subtract surrounding bins
[header_mX_snr, data_mX_snr] = RLW_SNR(header_mX, data_mX,...
            'xstart', snr_bins(1), 'xend', snr_bins(2)); 

header_mX_snr.name = sprintf('%s_snr-%d-%d_mX', ...
                             fname, snr_bins(1), snr_bins(2)); 
                         
if do_save
    CLW_save(save_path, header_mX_snr, data_mX_snr); 
end

% ---- complex FFT -----

[header_X, data_X] = RLW_FFT(header, data, 'output', 'complex'); 

[header_X, data_X] = crop_lw_data(header_X, data_X, fmax); 

data_X(:, :, :, :, :, 1) = 0; 

header_X.name = sprintf('%s_snr-none_X', fname);
            
if do_save
    CLW_save(save_path, header_X, data_X); 
end

% get noise-subtracted complex spectra 
header_X_snr = header_X; 

[data_X_snr] = subtract_noise_bins_complex(...
                    data_X, snr_bins(1), snr_bins(2)); 
               
header_X_snr.name = sprintf('%s_snr-%d-%d_X', ...
                        fname, snr_bins(1), snr_bins(2)); ; 

if do_save
    CLW_save(save_path, header_X_snr, data_X_snr); 
end


out = []; 

out.header_mX = header_mX; 
out.data_mX = data_mX; 

out.header_mX_snr = header_mX_snr; 
out.data_mX_snr = data_mX_snr; 

out.header_X = header_X; 
out.data_X = data_X; 

out.header_X_snr = header_X_snr; 
out.data_X_snr = data_X_snr; 











      
      
      
      