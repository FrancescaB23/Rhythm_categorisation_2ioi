function [header, data] = merge_datasets(datasets) 

header = datasets(1).header; 
data = datasets(1).data; 

for i_datset=2:length(datasets)
    
    for i_event=1:length(datasets(i_datset).header.events)
        
        header.events(end+1).latency = ...
            (header.datasize(end)*header.xstep) + ...
            datasets(i_datset).header.events(i_event).latency; 
        
        header.events(end).epoch = 1; 
        
        header.events(end).code = ...
            datasets(i_datset).header.events(i_event).code; 
    end
    
    data = cat(6,data,datasets(i_datset).data); 
    
    header.datasize(end) = header.datasize(end) + ...
                                datasets(i_datset).header.datasize(end); 

end
