function [Props, setup] = fx_read_setup(my_tdms_struct)
% Addme read the experiment setup of tdms. file from Cscan
%   my_tdms_struct: the my_tdms_struct read by TDMS_getStruct(filename);
Props = my_tdms_struct.Props;
Settings = Props.Settings;
newStr = split( Settings , '<Name>');
% find out the values
expression = '>\w*\.*\w*(\[\w*\])*<';
[values, ~] = regexp(newStr, expression, 'match', 'tokens');
% find out the names
expression = '(\w*\s)*\w*(\[.*\])*(\(.*\))*<';
[names, ~] = regexp(newStr, expression, 'match', 'tokens');
% make the dict 
setup = struct("Settings", 1);
for i = 2:length(names)
    filed = names{i}{1, 1};
    filed = erase(filed, {'<', ' ', ']', '[', '/', '(', ')', ':', '='});
    value = values{i, 1};
    value = erase(value, '<');
    value = erase(value, '>');
    setup.(filed) = value;
end

% save 
% save('Setups.mat', 'Props', 'setup');

end

