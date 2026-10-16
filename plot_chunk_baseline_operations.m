function [snr_f, snr_chunk_plot, zbl_f, zbl_chunk_plot, ...
    sbl_f, sbl_chunk_plot] = ...
    plot_chunk_baseline_operations(mX, freq, frex, bin_from, bin_to)

% Returns baseline operations on a chunk of frequencies centered around the
% frequency(ies) of interest 
%
% operations in the script:
%   1) Uncorrected amplitude frequency values are chunked around the
%   frequencies of interest
%   2) chunks are summed
%   3) baseline operations are performed: SNR, zSNR (z-scores), baseline
%   correction
%   These operations are performed on all frequencies included in the
%   [-bin_to, bin_to] range for visualisation purposes. Hence, a wider
%   chunk is at first considered in order to have enough data points to
%   perform the same operations also on bins that constitute the baseline
%   for the frequencies of interest.
%
%
% Inputs:
% mX : raw magnitude spectra with frequency as the last dimension (letswave
%       output). for the script to work, all other dimensions must be equal
%       to 1 (i.e., average across conditions and channels)
% freq : array_like
%     Frequencies for the FFT. 
% frex : array_like
%     Frequencies of interest [Hz]
% bin_from : int
%     Closest frequency bin (on both sides) used to calculate noise properties.
% bin_to : int
%     Furthest frequency bin (on both sides) used to calculate noise properties.
%
%
% Outputs:
% snr_f: snr at the frequency of interest
% snr_chunk_plot: snr values of all frequencies centered around frequency
%        of interest in the interval [-bin_to, bin_to]
% zbl_f: zbl (zSNR, zscore) at the frequency of interest
% zbl_chunk_plot: zbl values of all frequencies centered around frequency
%        of interest in the interval [-bin_to, bin_to]
% sbl_f: baseline corrected amplitude at the frequency of interest
% sbl_chunk_plot: baseline corrected amplitude values of all frequencies 
%        centered around frequency of interest in the interval 
%        [-bin_to, bin_to]


% if all but last dimension of mX are not equal to 1, error
size_mX = size(mX);
if ~all(size_mX(1:5) == 1)
    error('error: the first 5 dimensions of the data mX must be equal to 1');
end


% bins index of the frequencies of interest
frex_idx = dsearchn(freq', frex'); 


% extract chunk centered around frequencies of interest and sum them.
% Enough points should be included so to be able to perform baseline
% operations in the [-bin_to, bin_to] range

for i_f = 1:length(frex_idx)

    idx_start = frex_idx(i_f) - 2 * bin_to;
    idx_end = frex_idx(i_f) + 2 * bin_to;

    if idx_start <= 0
        error('Frex %d is too close to the start of the input (idx_start = %d).', ...
              i_f, idx_start); 
    end
    if idx_end > size(mX, 6)
        error('Frex %d is too close to the end of the input (idx_end = %d).', ...
              i_f, idx_end); 
    end

    chunk_per_f(i_f, :) = squeeze(mX(:, :, :, :, :, (idx_start : idx_end)))';
end

% sum data over chunks centered around frequencies of interest
chunk = sum(chunk_per_f, 1);


% calculate baseline operations at all bins in the [-bin_to, +bin_to] range

snr_chunk = zeros(1, (4*bin_to)+1);
zbl_chunk = zeros(1, (4*bin_to)+1);
sbl_chunk = zeros(1, (4*bin_to)+1);

for j = bin_to+1 : (bin_to*3)+1 % indices referring to frequencies in the [-bin_to, +bin_to] range

    % indices baseline bins for bin j
    idx_bl = [j-bin_to : j-bin_from j+bin_from : j+bin_to];

    % compute baseline and its standard deviation for bin j
    bl = mean(chunk(idx_bl));
    stdbl = std(chunk(idx_bl));

    snr_chunk(j) = chunk(j)/bl;
    zbl_chunk(j) = (chunk(j)-bl)/stdbl;
    sbl_chunk(j) = chunk(j)-bl;

end

% only keep bins relevant for the plot
snr_chunk_plot = snr_chunk(bin_to+1 : (bin_to*3)+1);
snr_f = snr_chunk_plot(bin_to+1);

zbl_chunk_plot = zbl_chunk(bin_to+1 : (bin_to*3)+1);
zbl_f = zbl_chunk_plot(bin_to+1);

sbl_chunk_plot = sbl_chunk(bin_to+1 : (bin_to*3)+1);
sbl_f = sbl_chunk_plot(bin_to+1);

end