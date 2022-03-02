function cf = fx_plot_estthick_tightplot(t_space, techs, errors, stds, usthickness, std_mean_threshold)
% plot 9 errors bar in 3*3 figure tightly
% snrs:     3 levels of snr
% techs:   3 kinds of techs
% errors: 3*3 cell of errors
% stds:    3*3 cell of stds
% usthickness: uniformly-spaced thickness
% plot in errorbars
cf                = figure('Name', ['best_results_' 'estimated_thickness']); % modify this 
set(cf, 'Position', [-200, -200, 1100, 700], 'color', 'white');  % modify this 
cas               = cell(length(techs), 1);
techs_num         = length(techs);
figlabel          = {'(a)', '(b)', '(c)'};
for j = 1:length(techs)
    cas{j}         = subplot(1, techs_num, j);
    ca2            = cas{j};
    mean_thickness = errors{j};
    std_thickness  = stds{j};
    untrackable    = mean_thickness>usthickness*(1+std_mean_threshold(1)) | ...
        mean_thickness<usthickness*(1-std_mean_threshold(1)) | ...
        std_thickness>usthickness*std_mean_threshold(2);
%     % the positions after it are set as untrackable
%     k                  = find(untrackable); 
%     untrackable(k:end) = 1;
    %
    h = fx_shaded_errorbars(ca2, untrackable, t_space, mean_thickness, std_thickness);
    hold on;
    yl = yline(usthickness, 'k--', '', 'LineWidth',2);
    yl.LabelHorizontalAlignment = 'right';
    grid on;
    % legend
    legend({'Estimated ply thickness', 'Uniformly-spaced ply thickness'}, 'location', 'best');
    %             xticks(1:2:24);
    xlabel('\fontname {times new roman} layer i', 'fontsize', 16);
    %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
    set(ca2, 'Fontname', 'times new Roman', 'Fontsize', 16);
    set(ca2, 'linewidth', 2);
    if j==1 % first col
        ylabel('Thickness (mm)', 'fontsize', 16, 'fontname', 'times new roman');
    elseif j==3 % last col
        set(ca2, 'YTickLabel', []);
%         yl = yline(usthickness, 'k--', ...
%             'uniformly-spaced thickness', 'LineWidth',2);
%         yl.LabelHorizontalAlignment = 'right';
    else
        set(ca2, 'YTickLabel',[]);
    end
    % adjust the pos.
    p       = get(ca2, 'Position');
    p(1)    = p(1) - 0.03;
    p(2)    = p(2) + 0.02;
    p(3)    = p(3) + 0.05;
%     p(4)    = p(4) + 0.07;
    set(ca2, 'Position', p);
    title(techs{j}, 'fontsize', 16, 'fontname', 'times new roman');
    xlabel({'\fontname{times new roman} Ply p', figlabel{j}}, 'fontsize', 16);
    xlim([t_space(1)-0.5, t_space(end)+0.5]);
    grid on;
end

for j = 1:length(cas)-1
    linkaxes([cas{j}, cas{j+1}], 'xy');
end

% ca               = subplot(3, 3, 1);
% errorbar(ca, t_space, Origin_errors_50*100, Origin_std_50*100, '-s',...
%     'MarkerSize', 5, 'Linewidth', 2,'MarkerEdgeColor','red','MarkerFaceColor','red');
% xlabel({'\fontname{times new roman} TOF (\mus)', '(a)'}, 'fontsize', 16);
% ylabel('\fontname{times new roman} Errors (mean and std)', 'fontsize', 16);
% xlim([t_space(1) t_space(end)]);
% yl                             = ylim;
% ytickformat('percentage');
% set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
% set(ca, 'linewidth', 2);
% grid on;
% % adjust pos.
% set(ca,'XTickLabel',[]);
% p         = get(ca, 'Position');
% p(3)    = p(3) + 0.05;
% set(ca, 'Position', p);
% ca1               = subplot(3, 3, 2);
% errorbar(ca1, t_space, WDARextrap_error*100, WDARextrap_std*100, '-s',...
%     'MarkerSize', 5, 'Linewidth', 2,'MarkerEdgeColor','red','MarkerFaceColor','red');
% xlabel({'\fontname{times new roman} TOF (\mus)', '(b)'}, 'fontsize', 16);
% %             ylabel('\fontname{times new roman} Errors (mean and std)', 'fontsize', 16);
% xlim([t_space(1) t_space(end)]);
% ylim(yl);
% ytickformat('percentage');
% set(ca1, 'Fontname', 'times new Roman', 'Fontsize', 16);
% set(ca1, 'linewidth', 2);
% grid on;
% set(ca1, 'YTickLabel', []);
% % adjust pos.
% set(ca1,'XTickLabel',[]);
% p         = get(ca1, 'Position');
% p(3)    = p(3) + 0.05;
% set(ca1, 'Position', p);
% ca2               = subplot(1, 3, 3);
% errorbar(ca2, t_space, AS_LG_errors*100, AS_LG_std*100, '-s',...
%     'MarkerSize', 5, 'Linewidth', 2,'MarkerEdgeColor','red','MarkerFaceColor','red');
% xlabel({'\fontname{times new roman} TOF (\mus)', '(c)'}, 'fontsize', 16);
% %             ylabel('\fontname{times new roman} Errors (mean and std)', 'fontsize', 16);
% xlim([t_space(1) t_space(end)]);
% ylim(yl);
% ytickformat('percentage');
% set(ca2, 'Fontname', 'times new Roman', 'Fontsize', 16);
% set(ca2, 'linewidth', 2);
% grid on;
% set(ca2, 'YTickLabel', []);
% % adjust pos.
% set(ca2,'XTickLabel',[]);
% p         = get(ca2, 'Position');
% p(3)    = p(3) + 0.05;
% set(ca2, 'Position', p);
% hold on;
% linkaxes([ca ca1 ca2],'xy');
            
end

