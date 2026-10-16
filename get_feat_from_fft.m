function [feat_all_cond] = get_feat_from_fft(header, data, freq, frex)

% compute feature vector from magnitude or complex fft (specify)
% the output is a matrix where each column is a feature vector for one
% condition. Adapted from get_feat_from_X to use as input magnitude or
% complex FFT
%
% Parameters:
%
% header:   header, letswave structure
% data:     data, letswave structure, magnitude or complex fft
% freq:     array, vector containing the frequencies at each time point
% frex:     array, frequencies of interest in Hz


%%

% sanity check that the input is actually complex/amplitude spectrum as
% specified in the header

switch header.filetype

    case 'frequency_complex'

        if isreal(data)
            error('seems like you provided magnitude spectrum, but we need complex');
        end

    case 'frequency_amplitude'

        if ~isreal(data)
            error('seems like you provided complex spectrum, but we need magnitude');
        end

    otherwise
        error('unexpected header.filetype: %s', header.filetype);

end


%%
% indices of frequencies of interest to extract data for feature vector
frex_idx = dsearchn(ensure_col(freq), ensure_col(frex))'; 

n_cond = size(data, 1);
n_chans = size(data, 2);


switch header.filetype

    case 'frequency_complex'


        for i_cond = 1:n_cond

            feat_cond = [];

            for i_freq = 1:length(frex_idx)

                feat_cond = [feat_cond; real(squeeze(data(i_cond, 1:n_chans, 1, 1, 1, frex_idx(i_freq))))'];
            end

            for i_freq = 1:length(frex_idx)

                feat_cond = [feat_cond; imag(squeeze(data(i_cond, 1:n_chans, 1, 1, 1, frex_idx(i_freq))))'];
            end

            feat_all_cond(:, i_cond) = feat_cond;

        end

        expected_size = [2*length(frex_idx)*n_chans, n_cond];


    case 'frequency_amplitude'

        for i_cond = 1:n_cond

            feat_cond = [];

            for i_freq = 1:length(frex_idx)

                feat_cond = [feat_cond; squeeze(data(i_cond, 1:n_chans, 1, 1, 1, frex_idx(i_freq)))'];
            end

            feat_all_cond(:, i_cond) = feat_cond;

        end

        expected_size = [length(frex_idx)*n_chans, n_cond];
end



if ~isequal(size(feat_all_cond), expected_size)

    error('feat_all_cond is not of the expected size')
end


end