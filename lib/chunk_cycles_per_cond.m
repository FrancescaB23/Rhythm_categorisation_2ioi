function [header, data] = chunk_cycles_per_cond(header_trial, data_trial, cycle_dur, experiment)

n_cond = length(unique({header_trial.events.code})); 
    
events = []; 
data = []; 
c_ep = 1; 

data_all = {}; 

for i_cond=1:n_cond

    % get only epochs for this condition 
    if strcmp(experiment, 'newborn')
        ep_mask = strcmp({header_trial.events.code}, sprintf('S%03d', i_cond));
    else
        ep_mask = strcmp({header_trial.events.code}, num2str(i_cond)); 
    end

    [header_cond, data_cond] = RLW_arrange_epochs(...
                                header_trial, data_trial, find(ep_mask)); 

    data_cond_chunk = epoch_chunks(data_cond, 1/header_cond.xstep, cycle_dur); 
    
    for i_trial=1:size(data_cond_chunk, 2)
        
        for i_chunk=1:size(data_cond_chunk, 1)
        
            new_event = struct('code', num2str(i_cond), ...
                               'latency', 0, ...
                               'epoch', c_ep, ...
                               'trial', i_trial, ...
                               'cycle', i_chunk); 

            events = [events, new_event]; 

            c_ep = c_ep + 1; 
            
        end
    end
    
    data_all = [data_all, ...
        reshape(data_cond_chunk, ...
                [], header_cond.datasize(2), 1, 1, 1, size(data_cond_chunk, 7)...
                )...
                ]; 
    
end

data = cat(1, data_all{:}); 

header = header_trial;
header.datasize = size(data); 
header.events = events; 



% % sanity check 
% i_cond = 13; 
% i_trial = 6; 
% i_cycle = 30; 
% 
% ep_mask = strcmp({header_trial.events.code}, num2str(i_cond)); 
% [header_cond, data_cond] = RLW_arrange_epochs(header_trial, data_trial, find(ep_mask)); 
% data_cond_chunk = epoch_chunks(data_cond, 1/header_cond.xstep, par.cycle_dur); 
% x1 = squeeze(data_cond_chunk(i_cycle, i_trial, 1, 1, 1, 1, :)); 
% 
% chunk_mask = strcmp({events.code}, num2str(i_cond)) &...
%     [events.trial] == i_trial &...
%     [events.cycle] == i_cycle; 
% 
% x2 = squeeze(data(find(chunk_mask), 1, 1, 1, 1, :)); 
% 
% figure 
% plot(x1, 'linew', 2)
% hold on 
% plot(x2, 'linew', 2)
% 
% assert(isequal(x1, x2))

