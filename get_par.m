function par = get_par(varargin)

parser = inputParser(); 

addParameter(parser, 'experiment', '2ioi_conga'); % 2ioi_conga, newborn, 2ioi_stem, 3ioi_nonsweep

parse(parser, varargin{:}); 

experiment = parser.Results.experiment; 

%% prepare paths, add libraries

[~, hostname] = system('hostname');
hostname = deblank(hostname); 

if strcmp(hostname, 'mac-BX23-117.local')
    
    % rnb tools 
    rnb_tools_path = '~/projects_git/rnb_tools/src'; 
    
    % ftrsa
    ftrsa_path = '/Users/tomaslenc/projects_git/ftrsa/src'; 

    % letswave 6
    lw6_path = '~/projects_git/letswave6'; 

    % letswave 
    lw7_path = '~/projects_git/letswave7'; 

    % path to root of the data folder 
    if strcmp(experiment, '2ioi_stem')
        data_path = '~/projects_backed_up/Rhytotype/data_eeg_2ioi_stem'; 
    elseif strcmp(experiment, '2ioi')
        data_path = '~/projects_backed_up/Rhytotype/data_eeg_2ioi'; 
    elseif strcmp(experiment, '3ioi_nonsweep')
        data_path = '~/projects_backed_up/Rhytotype/data_3ioi_nonsweep'; 
    else
        error('experiment "%s" not recognized...', experiment); 
    end
   
elseif contains(hostname, 'tomo-office-desktop')
    
    % rnb tools 
    rnb_tools_path = '~/Documents/rnb_tools/src'; 
    
    % ftrsa
    ftrsa_path = '~/Documents/ftrsa/src'; 

    % letswave 6
    lw6_path = '~/Documents/letswave6'; 

    % letswave 
    lw7_path = '~/Documents/letswave7'; 

    % path to root of the data folder 
    if strcmp(experiment, '2ioi_stem')
        data_path = '/DATA1/Rhototype_analysis/data_2ioi_stem'; 
    elseif strcmp(experiment, '2ioi')
        data_path = '/DATA1/Rhototype_analysis/data_2ioi'; 
    elseif strcmp(experiment, '3ioi_nonsweep')
        data_path = ''; 
    else
        error('experiment "%s" not recognized...', experiment); 
    end

elseif contains(hostname, 'tux')
    
    % rnb tools 
    rnb_tools_path = '/datadisk/projects_git_dl/rnb_tools/src'; 

    % ftrsa
    ftrsa_path = ''; 

    % letswave 6
    lw6_path = '/home/tomo/Documents/MATLAB/letswave6'; 

    % letswave 
    lw7_path = '/home/tomo/Documents/MATLAB/letswave7'; 

    % path to root of the data folder 
    if strcmp(experiment, '2ioi_stem')
        data_path = '/datadisk/projects_backed_up/Rhytotype/data_eeg_2ioi'; 
    elseif strcmp(experiment, '2ioi')
        data_path = '/datadisk/projects_backed_up/Rhytotype/data_eeg'; 
    else
        error('experiment "%s" not recognized...', experiment); 
    end  

elseif strcmp(hostname, 'mac-BX22-227.local')

    % rnb tools 
    rnb_tools_path = ['/Users/barbero/Documents/GitHub/' ...
        'rnb_tools_fblaptop/src']; 
    
    % ftrsa
    ftrsa_path = '/Users/barbero/Documents/GitHub/ftrsa/src'; 

    % artefact blocking algorithm
    ab_path = ['/Users/barbero/Documents/GitHub/' ...
        'Rhototype_2ioi_newborns/library/AB'];

    % letswave 6
    lw6_path = ['/Users/barbero/Documents/GitHub/' ...
        'Rhototype_2ioi_newborns/library/letswave6-master_custom'];

    % letswave 
    lw7_path = ['/Users/barbero/Documents/GitHub/' ...
        'Rhototype_2ioi_newborns/library/letswave7-master'];

    % path to root of the data folder 
    if strcmp(experiment, 'newborn')
        data_path = ['/Users/barbero/Documents/Projects/' ...
            'Rhototype_2ioi_newborns/newborns/data'];
    elseif strcmp(experiment, '2ioi_conga')
        data_path = ['/Users/barbero/Documents/Projects/' ...
            'Rhototype_2ioi_newborns/adults/data'];
    else
        error('experiment "%s" not recognized...', experiment); 
    end
    
    
end

addpath(genpath(rnb_tools_path)); 
addpath(genpath(ftrsa_path)); 
addpath(genpath(lw6_path)); 
addpath(genpath('lib')); 

% paths to data subfolders 
raw_path = fullfile(data_path, 'raw'); 
deriv_path = fullfile(data_path, 'derivatives'); 
results_path = fullfile(deriv_path, 'results'); 
summary_path = fullfile(deriv_path, 'summary'); 

if ~isdir(raw_path)
    mkdir(raw_path)
end
if ~isdir(deriv_path)
    mkdir(deriv_path)
end
if ~isdir(summary_path)
    mkdir(summary_path)
end

%% experiment-specific parameters

if strcmp(experiment, 'newborn')

    % trial duration 
    trial_dur = 30;
    % rhythmic pattern duration 
    cycle_dur = 0.750; 
    % number of conditions 
    n_cond = 10;

elseif strcmp(experiment, '2ioi_conga')

    trial_dur = 30;
    cycle_dur = 0.750; 
    n_cond = 10;

    % parameters to retrieve tapping force from six force signal of the
    % tapping box
    tap_baseline = 1;
    % calibration coefficients tapping box, Volt to force. Calibration
    % matrix from December 2023
    calib_coeff_Fz = [10.52079  -0.21329  10.53074 -0.22406  10.62495 ...
                    -0.31096];

else
    error('experiment "%s" not recognized...', experiment); 
end


%%

% time left before start and after the end of each trial (to allow later
% filtering without edge artifacts)
epoch_buffer_dur = 5; 

% number of channels to use for interpolation 
n_chans_interp = 3; 

if strcmp(experiment, 'newborn')
    filt_cut_low = 0.5;
    filt_cut_high = 30; 
else
    filt_cut_low = 0.1; 
    filt_cut_high = 64; 
end
filt_order_low = 4; 
filt_order_high = 4; 

filt_ica_cut_low = 1; 
filt_ica_cut_high = 64; 
filt_ica_order_low = 4; 
filt_ica_order_high = 4; 

decim_factor = 4; 
decim_factor_ica = 4; 

filt_erp_cut_low = 10; 
filt_tap_erp_cut_low = 10; 


%% Channel parameters

if strcmp(experiment, 'newborn')

    % Electrodes of the outer circle (to discard)
    chans_outer_circle = {'E17' 'E43' 'E48' 'E49' 'E56' 'E63' 'E68' ...
        'E73' 'E81' 'E88' 'E94' 'E99' 'E107' 'E113' 'E119' 'E120' ...
        'E125' 'E126' 'E127' 'E128'};

    eeg_ref_chans = {'E57', 'E100'}; 

%     % Channels labels - excluding three outer circles including mastoids
%     chans_internal = {'E3' 'E4' 'E5' 'E6' 'E7' 'E10' 'E11' 'E12' ...
%         'E13' 'E16' 'E18' 'E19' 'E20' 'E23' 'E24' 'E27' 'E28' 'E29' ...
%         'E30' 'E31' 'E34' 'E35' 'E36' 'E37' 'E40' 'E41' 'E42' 'E46' ...
%         'E47' 'E51' 'E52' 'E53' 'E54' 'E55' 'E59' 'E60' 'E61' 'E62' ...
%         'E66' 'E67' 'E71' 'E72' 'E76' 'E77' 'E78' 'E79' 'E80' 'E84' ...
%         'E85' 'E86' 'E87' 'E91' 'E92' 'E93' 'E97' 'E98' 'E102' 'E103' ...
%         'E104' 'E105' 'E106' 'E109' 'E110' 'E111' 'E112' 'E116' 'E117' ...
%         'E118' 'E123' 'E124'};

else

%     roi_front_chans = {
%             'F1'
%             'Fz'
%             'F2'
%             'FC1'
%             'FCz'
%             'FC2'
%             'C1'
%             'Cz'
%             'C2'
%         }; 
    
    eeg_ref_chans = {'TP9', 'TP10'}; 

end


%%

% number of frequency bins for SNR subtraction and zSNR estimation 
snr_bins = [2, 5]; 

% number of permutations when testing fit of a categorical model at
% individual participant level 
if strcmp(experiment, '2ioi')
    n_perm_rsa_individual = 5000; 
else
    n_perm_rsa_individual = 10000; 
end

% a significant category has to have more pixels (i.e. conditions) than
% this
n_mini_categ_skip = 1; 

%% return structure 

clear varargin

w = whos;
par = []; 
for a = 1:length(w) 
    par.(w(a).name) = eval(w(a).name); 
end