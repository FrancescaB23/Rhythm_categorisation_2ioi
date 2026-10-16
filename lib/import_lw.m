function letswaveNewPath = import_lw(varargin)

par = get_par(); 

% default (or 6): import LW6
% pass 7 to import LW7
% pass -1 to remove everything from path
lw6 = par.lw6_path;  
lw7 = par.lw7_path; 

if nargin==1 && varargin{1}==-1
    warning('off')
    rmpath(genpath(lw6)); 
    rmpath(genpath(lw7)); 
    fprintf('\nremoving everything from path...\n\n'); 
    warning('on')
    return
end

if nargin==1 && varargin{1}==7
    letswaveOldPath = lw6; 
    letswaveNewPath = lw7; 
    fprintf('\nadding LW7 to path...\n\n'); 
else
    letswaveOldPath = lw7; 
    letswaveNewPath = lw6; 
    fprintf('\nadding LW6 to path...\n\n'); 
end

warning('off');
rmpath(genpath(letswaveOldPath));
addpath(genpath(letswaveNewPath));
warning('on');

