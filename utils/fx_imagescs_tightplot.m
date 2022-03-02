function cf = fx_imagescs_tightplot(X, Y, layers, techs, profile_layer)
% plot 9 images in 3*3 figure tightly
% X, Y: the x, y axis ticks
% layers:  3 kinds of layers, cell
% techs:   3 kinds of techs
% profile_layer: 3*3 cell of images

% plot in errorbars
cf               = figure('Name', ['best_results_' 'profiles']);
set(cf, 'Position', [-300, -300, 1200, 1200], 'color', 'white');
cas              = cell(9);
for i = 1:3
    for j = 1:3
        cas{i*3+j-3}           = subplot(3, 3, i*3+j-3);
        imagesc(cas{i*3+j-3}, X, Y, profile_layer{i, j});
        colormap('jet');
        hold on;
        set(cas{i*3+j-3}, 'Fontname', 'times new Roman', 'Fontsize', 16);
        set(cas{i*3+j-3}, 'linewidth', 2);
        if j==1 % first col.
            ylabel({layers{i}, 'Y displacement (mm)'}, 'fontsize', 16, 'fontname', 'times new roman');
            cl = caxis;
            % adjust the pos.
            p       = get(cas{i*3+j-3}, 'Position');
            p(2)    = p(2) - 0.03;
            p(3)    = p(3) + 0.04;  
            p(4)    = p(4) + 0.07;
            set(cas{i*3+j-3}, 'Position', p);
        elseif j==3 % last col.
            caxis(cl);
            h      = colorbar;
            pcb    = h.Position; %gets the positon and size of the color bar
            pcb(4) = pcb(4);
            set(h,'Position', pcb)% To change size
            if i==1 % last col. and first row
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Depth (mm)');
            end
            set(cas{i*3+j-3}, 'YTickLabel', []);
            p       = get(cas{i*3+j-3}, 'Position');
            p(1)    = p(1) - 0.04;
            p(2)    = p(2) - 0.03;
            p(3)    = p(3) + 0.07;
            p(4)    = p(4) + 0.07;
            set(cas{i*3+j-3}, 'Position', p);
        else
            caxis(cl);
            set(cas{i*3+j-3}, 'YTickLabel', []);
            % adjust the pos.
            p       = get(cas{i*3+j-3}, 'Position');
            p(1)    = p(1) - 0.02;
            p(2)    = p(2) - 0.03;
            p(3)    = p(3) + 0.04;  
            p(4)    = p(4) + 0.07;
            set(cas{i*3+j-3}, 'Position', p);
        end
        if i==1 % first row 
            title(techs{j}, 'fontsize', 16, 'fontname', 'times new roman');
            set(cas{i*3+j-3}, 'XTickLabel', []);
            caxis([1 1.3]);
        elseif i==3 % last row
            xlabel({'\fontname{times new roman} X displacement (mm) '}, 'fontsize', 16);
            caxis([4.9 5.2]);
        else
            set(cas{i*3+j-3}, 'XTickLabel', []);
            caxis([3.3 3.6]);
        end
    end
end

% linkaxes([cas{1}, cas{2}, cas{3}, ...
%     cas{4}, cas{5}, cas{6}, ...
%     cas{7}, cas{8}, cas{9}], 'xy');

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

