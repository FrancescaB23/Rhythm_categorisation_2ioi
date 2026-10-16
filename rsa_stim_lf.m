function rdm = rsa_stim_lf(par, frex, varargin)

parser = inputParser(); 

addParameter(parser, 'method', 'pearson'); 
addParameter(parser, 'feat', 'X'); % X, mX, ratio

parse(parser, varargin{:}); 

method = parser.Results.method; 
feat_name = parser.Results.feat; 

%%

% load stimulus info 
load(fullfile(par.deriv_path, 'stimuli', 'stim.mat')); 

if strcmp(feat_name, 'ratio')
    
    ratios = stim.ratios(:, 1); 
    
    rdm = get_rdm_from_feat(ratios', 'method', 'L1'); 
    
else
    % get RDM from stimulus envelope 
    stim_env = cellfun(@(x) abs(hilbert(x)), stim.s, 'uni', 0); 

    [X_stim, freq_stim] = get_X(stim_env, stim.fs); 

    feat_stim = get_feat_from_X(X_stim, freq_stim, frex, 'feat', feat_name); 

    rdm = get_rdm_from_feat(feat_stim, 'method', method); 

end
