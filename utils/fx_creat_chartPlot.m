function fx_creat_chartPlot(maps, cData_collection, Stacking_sequence)
%%
% Number of color levels to create

nLevels = size(cData_collection, 2);
cf      = figure('Name', ['stacked_graymap_angleDistribute']);
set(cf, 'Position', [0, 0, 650, 1000], 'color', 'white');
% X data points for color patches
xData = [linspace(0, nLevels-1, nLevels); linspace(1, nLevels, nLevels); ...
    linspace(1, nLevels, nLevels); linspace(0, nLevels-1, nLevels)];

offset_para = 12;
% Create each color bar
% for uniforming the color scale 
c1 = max(cData_collection, [], 'all');
c2 = min(cData_collection, [], 'all');

for iMap = 1:length(maps)
    offset = offset_para*(length(maps) - iMap);
    yData = [zeros(2, nLevels); (offset_para-2)*ones(2, nLevels)] + offset;
    % Construct appropriate colormap.
    %     cData = cData_collection(iMap, :) * 100;
    % *** log ***
    %     cData = log10(cData_collection(iMap, :));
    % *** uniform ***
    cData = cData_collection(iMap, :) / max(cData_collection(iMap, :));
    % Display colormap chart
    patch('XData', xData, 'YData', yData, ...
        'EdgeColor', 'none', ...
        'FaceColor', 'flat', ...
        'FaceVertexCData', cData');
    colormap gray;
    %     set(gca,'ColorScale','log');
    % set limits for the caxis
    % *** log ***
    %     caxis([log10(c2) log10(c1)]);
    %     caxis([0 maxCdata]);
    rectangle(...
        'Position', [0-0.5, offset, nLevels+1, offset_para-2], ...
        'Curvature', [0 0], ...
        'EdgeColor', 'w',...
        'LineWidth', 1);
    text(0-65, offset, sprintf(' %s', maps{iMap}), ...
        'VerticalAlignment', 'bottom', ...
        'FontSize', 12, ...
        'FontName', 'Times New Romans');
    set(gca, 'Fontname', 'times new Roman', 'FontSize', 12);
    set(gca, 'linewidth', 2);
    line([0 0]+22 + 90 + Stacking_sequence(iMap), ...
        [offset offset+offset_para-2], 'color', 'red', 'linewidth', 1);
end
axis equal off
% title('Built-in Colormaps')
c = colorbar;
% set(gca,'ColorScale','log');
set(get(c, 'Title'), 'string', '\fontname {times new roman}\fontsize {12} Normalized value');           
% c.Ruler.TickLabelFormat='%g%%';
c.Location = 'northoutside';
c.FontSize = 12;
c.FontName = 'Times new romans';

line([0 0]+22, [-2 0], 'linewidth', 1, 'color', 'red');
line([0 0]+22+45, [-2 0], 'linewidth', 1, 'color', 'red');
line([0 0]+22+90, [-2 0], 'linewidth', 1, 'color', 'red');
line([0 0]+22+135, [-2 0], 'linewidth', 1, 'color', 'red');

end

