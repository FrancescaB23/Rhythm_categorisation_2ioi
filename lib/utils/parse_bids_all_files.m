function p = parse_bids_all_files(fpath, varargin)
% parse all filesnames in the given folder
parser = inputParser; 

addParameter(parser, 'incl_dirs', false); 
addParameter(parser, 'glob', '*'); 

parse(parser, varargin{:}); 

incl_dirs = parser.Results.incl_dirs; 
glob_str = parser.Results.glob; 

%%

d = dir(fullfile(fpath, glob_str)); 

if ~incl_dirs
    d = d(~[d.isdir]); 
end

p = cellfun(@(x) parse_bids_entities(x), {d.name}, 'uni', 0); 
p = cat(1, p{:}); 