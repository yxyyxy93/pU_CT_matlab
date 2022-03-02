function [fig, ax]= fx_showCS(fignum, inam_C_scan, z_index, name_colormap, titlestr)
% show b scan

switch nargin
    case 4
        titlestr = 'C scan section of '; 
    case 5
        disp(titlestr);    
end

% save and output the handles
fig = figure(fignum);

ax = imagesc(inam_C_scan);
hold on;

% figure setup
title([titlestr, num2str(z_index)], 'fontsize', 16, 'Fontname', 'times new Roman');   
colormap(name_colormap);
h = colorbar;
set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16}');
set(gca, 'Fontname', 'times new Roman');
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);
xlabel('\fontname {times new roman}', 'fontsize', 16);
ylabel('\fontname {times new roman}', 'fontsize', 16);
% draw the line indicating the front surface

set(gcf, 'Position', [0, 0, 800, 600], 'color', 'white');

end

