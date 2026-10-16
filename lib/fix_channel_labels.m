function [header, data] = fix_channel_labels(...
                        header, data, labels_path, lw6_path, varargin)

parser = inputParser; 

addParameter(parser, 'keep_eeg_only', true); 

parse(parser, varargin{:}); 

keep_eeg_only = parser.Results.keep_eeg_only; 

%% 

% remap channel labels to 10/20 
chan_labels = readtable(labels_path); 
chanlocs = header.chanlocs; 
for i_chan=1:length(chanlocs)
    old_lab = chanlocs(i_chan).labels; 
    idx = find(strcmpi(old_lab, chan_labels.old)); 
    if length(idx) == 1
        chanlocs(i_chan).labels = chan_labels.new{idx}; 
    end
end
header.chanlocs = chanlocs; 

if keep_eeg_only
    % only keep EEG channels
    [header, data] = RLW_arrange_channels(header, data, chan_labels.new); 
end

% edit electrode coordinates
loc_file_path = fullfile(lw6_path, ...
                         'core_functions', ...
                         'Standard-10-20-Cap81.locs');

locs = readlocs(loc_file_path);

for i_chan=1:length(header.chanlocs)

    idx_chan = find(strcmpi(header.chanlocs(i_chan).labels, {locs.labels}));

    if ~isempty(idx_chan)
        header.chanlocs(i_chan).theta           = locs(idx_chan).theta;
        header.chanlocs(i_chan).radius          = locs(idx_chan).radius;
        header.chanlocs(i_chan).sph_theta       = locs(idx_chan).sph_theta;
        header.chanlocs(i_chan).sph_phi         = locs(idx_chan).sph_phi;
        header.chanlocs(i_chan).sph_theta_besa  = locs(idx_chan).sph_theta_besa;
        header.chanlocs(i_chan).sph_phi_besa    = locs(idx_chan).sph_phi_besa;
        header.chanlocs(i_chan).X               = locs(idx_chan).X;
        header.chanlocs(i_chan).Y               = locs(idx_chan).Y;
        header.chanlocs(i_chan).Z               = locs(idx_chan).Z;
        header.chanlocs(i_chan).topo_enabled    = 1;
        header.chanlocs(i_chan).SEEG_enabled    = 0;
    end
end
