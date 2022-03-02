function  fx_plot_Bscan(ca, Exp_setup, AS_analysis, B_scan_phase, location_y, titlename)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here

set(ca, 'Fontname', 'times new Roman', 'fontsize', 16);
set(ca, 'linewidth', 1.5);
imagesc(ca, (1:Exp_setup.B_scan_points), AS_analysis.t * 1e6, B_scan_phase); 
% draw lines of the location of the interface
hold on;
colormap(ca, parula);
xlabel(ca, '\fontname {times new roman} scan index', 'fontsize', 16);
ylabel(ca, '\fontname {times new roman} time (\mus)', 'fontsize', 16);

ylim(ca, [-max(max(location_y)) * 0.1* 1e6 max(max(location_y)) *1.1* 1e6]);
h = colorbar(ca);
set(get(h,'Title'), 'string','\fontname {times new roman}\fontsize {16}');
title(ca, titlename, 'fontsize', 16, 'Fontname', 'times new Roman');

ca.FontSize = 12; 

end

