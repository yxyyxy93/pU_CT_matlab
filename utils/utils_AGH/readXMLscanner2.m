function scanSet = readXMLscanner2(fn);

%% Initialize variables.
filename = fn;
delimiter = '';

%% Format string for each line of text:
%   column1: text (%s)
% For more information, see the TEXTSCAN documentation.
formatSpec = '%s%[^\n\r]';

%% Open the text file.
fileID = fopen(filename,'r');

%% Read columns of data according to format string.
% This call is based on the structure of the file used to generate this
% code. If an error occurs for a different file, try regenerating the code
% from the Import Tool.
dataArray = textscan(fileID, formatSpec, 'Delimiter', delimiter, 'MultipleDelimsAsOne', true,  'ReturnOnError', false);

%% Close the text file.
fclose(fileID);

%% Post processing for unimportable data.
% No unimportable data rules were applied during the import, so no post
% processing code is included. To generate code which works for
% unimportable data, select unimportable cells in a file and regenerate the
% script.

%% Allocate imported array to column variable names
VarName2 = dataArray{:, 1};
%%
% for i =1:length(VarName2)
% disp([ num2str(i) ' ' VarName2{i}])
% end
selDat = [24 28 32 71];

ind1 = 1; clear scanerSet;
for j = 1 :length(selDat)
    napis = VarName2{selDat(j)}(7:(strfind(VarName2{selDat(j)},'</Name>')-1));
    wartosc = str2num(VarName2{selDat(j)+1}(6:(strfind(VarName2{selDat(j)+1},'</Val>')-1)));
    %     disp([ napis ' ' num2str(wartosc)])
    scanerSet(ind1).name = napis;
    scanerSet(ind1).val = wartosc;
    ind1 = ind1+1;
end
scanerSet(ind1).name = 'x range [mm]';
scanerSet(ind1).val = [str2num(VarName2{(45)+1}(6:(strfind(VarName2{(45)+1},'</Val>')-1)))...
    str2num(VarName2{(49)+1}(6:(strfind(VarName2{(49)+1},'</Val>')-1)))];
scanerSet(ind1+1).name = 'y range [mm]';
scanerSet(ind1+1).val = [str2num(VarName2{(59)+1}(6:(strfind(VarName2{(59)+1},'</Val>')-1)))...
    str2num(VarName2{(63)+1}(6:(strfind(VarName2{(63)+1},'</Val>')-1)))];


scanSet.yRange = abs(scanerSet(6).val(2)- scanerSet(6).val(1));     scanSet.xRange = abs(scanerSet(5).val(2)- scanerSet(5).val(1));
scanSet.yStep  = scanerSet(3).val;                                  scanSet.xSpeed = scanerSet(1).val;
