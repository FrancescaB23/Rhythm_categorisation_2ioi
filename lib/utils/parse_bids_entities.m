function p = parse_bids_entities(fname)  
% This is a parser for BIDS-like filename strcture that I adopted in my
% experiments. It's not valid BIDS!!! For example I often use several
% dashes and several characters in the entity values, which is not allowed
% in the BIDS specification. But life is too short to follow all rules...

[~, basename, ext] = fileparts(fname); 

p.entities = struct();
p.suffix = '';

% -Identify all the BIDS entity-label pairs present in the filename (delimited by "_")
[parts, dummy] = regexp(basename, '(?:_)+', 'split', 'match');

% Separate the entity from the label for each pair identified above
for i = 1:numel(parts)

%   [d, dummy] = regexp(parts{i}, '(?:\-)+', 'split', 'once', 'match');

  d = find(parts{i} == '-'); 

  switch length(d)

    case 0 % no - in entity, may be suffix
      p.suffix = parts{i};

    otherwise % normal entity
      p.entities.(parts{i}(1:d(1)-1)) = parts{i}(d(1)+1:end);

  end

end

