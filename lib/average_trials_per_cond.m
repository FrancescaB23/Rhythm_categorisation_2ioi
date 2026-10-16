function [header, data] = average_trials_per_cond(header_ep, data_ep, par, varargin)

parser = inputParser; 

addParameter(parser, 'polarity', 'sum'); % add, subtract
addParameter(parser, 'concatenate_epochs', false); % add, subtract
addParameter(parser, 'par', []); 

parse(parser, varargin{:}); 

polarity = parser.Results.polarity; 
concatenate_epochs = parser.Results.concatenate_epochs; 


% if parameter structure hasn't been explicitly passed, load a new one from
% the file - this allows the possibility for the user to pass an adjusted
% par structure intead of the default one
if isempty(par)
    par = get_par(); 
end


%%

n_cond = length(unique({header_ep.events.code})); 
      
datasets = []; 

for i_cond=1:n_cond

    % get only epochs for this condition 
    if strcmp(par.experiment, '2ioi_conga')
        ep_mask = strcmp({header_ep.events.code}, num2str(i_cond)); 
    elseif strcmp(par.experiment, 'newborn')
        ep_mask = strcmp({header_ep.events.code}, sprintf('S0%02d',i_cond)); % for newborn data
    end

    [header_cond, data_cond] = RLW_arrange_epochs(header_ep, data_ep, find(ep_mask)); 

    % averrge trials 
    if concatenate_epochs

        [header_cond, data_cond] = RLW_concatenate_epochs(...
                        header_cond, data_cond, [1:header_cond.datasize(1)]); 

        warning('concatenating epochs');

    else

        if strcmp(polarity, 'subtract')
            
            warning('subtracting positive and negative polarity trials'); 
            
            if ~isfield(header_cond.events, 'polarity')
                error('no polarity information in the header'); 
            end
            
            polarity_per_epoch = [header_cond.events.polarity]; 
            
            for i_ep=1:size(data_cond, 1)
                data_cond(i_ep, :, 1, 1, 1, :) = data_cond(i_ep, :, 1, 1, 1, :) ...
                                                    * polarity_per_epoch(i_ep); 
            end
            
        end

        [header_cond, data_cond] = RLW_average_epochs(header_cond, data_cond); 

        header_cond.events = header_cond.events(1); 
        header_cond.events(1).epoch = 1; 
    end
    
    datasets(i_cond).header = header_cond; 
    datasets(i_cond).data = data_cond; 

end

%%

[header, data] = RLW_merge_epochs(datasets, [1:n_cond]); 



