function h = fx_shaded_errorbars(ca, untrackable, t_space, Origin_mean, Origin_std)

% plot the errorbars by shading the untrackable points
alpha   = 0.3;
h       = errorbar(ca, t_space(untrackable), Origin_mean(untrackable), ...
    Origin_std(untrackable), 's', 'MarkerSize',5, 'MarkerFaceColor', [1 1 1]*(1-alpha), ...
    'CapSize', 10, 'LineWidth', 2, 'LineStyle', 'none', 'HandleVisibility','off'); % no showing on legend
h.Color = [1 1 1]*(1-alpha);
% Set transparency
set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', ...
    'ColorData', [h.Line.ColorData(1:3); 255*alpha]);
hold on;
h       = errorbar(ca, t_space(~untrackable), Origin_mean(~untrackable), ...
    Origin_std(~untrackable), 's', 'MarkerSize', 5, 'MarkerEdgeColor', 'red', ...
    'MarkerFaceColor','red', 'CapSize', 10, 'LineWidth', 2, 'LineStyle', 'none');
h.Color = [0.00,0.45,0.74];

end

