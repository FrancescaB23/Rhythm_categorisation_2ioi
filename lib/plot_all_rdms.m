function plot_all_rdms(pnl, subjects, rdms, best_models, pvals, z_snrs)

pnl.pack('v', length(subjects));     

for i_sub=1:length(subjects)
            
    pnl(i_sub).pack('h', 2); 
       
    ax = pnl(i_sub, 1).select(); 
    plot_rdm(rdms{i_sub}, 'ax', ax); 
    ax.YDir = 'reverse'; 
    
    ax = pnl(i_sub, 2).select(); 
    plot_rdm(best_models{i_sub}, 'ax', ax); 
    ax.YDir = 'reverse'; 
    
    tit = sprintf('p=%.1g', pvals(i_sub)); 
    if nargin == 6
       tit = sprintf('%s, z-snr=%.2f', tit, z_snrs(i_sub)); 
    end
    pnl(i_sub).title(tit); 
    
    pnl(i_sub).ylabel(sprintf('%03d', subjects(i_sub))); 

end

pnl.de.margin = [0, 0, 0, 8]; 
pnl.marginright = 10; 

