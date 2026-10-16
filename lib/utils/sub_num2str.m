function sub_str=sub_num2str(sub)

if isnumeric(sub)
    sub_str = sprintf('sub-%03d', sub); 
elseif ischar(sub)
    sub_str = sprintf('sub-%s', sub); 
end