function cf = fx_plot_bscans_tightplot(xs, zs, Bscans_bg, Bscans_track, front_line, rear_line, techs, fss, fxs, usTOF)
% plot 9 errors bar in 3*3 figure tightly
% X_space, Y_space: x, y axises
% techs:   3 kinds of techs
% Bscans_bg:  the 2d bscan images of background
% Bscans_track:  the 2d bscan images of interply tracks
% front_line, rear_line: indexes of the front and back surface;
% fss, fxs: sampling rate of z and x
% plot in errorbars
cf                = figure('Name', ['best_results_' 'Bscans_interplies']); % modify this 
set(cf, 'Position', [-200, -200, 1400, 700], 'color', 'white');  % modify this 
cas               = cell(length(techs), 1);
techs_num         = length(techs);
figlabel          = {'(a)', '(b)', '(c)'};
for j = 1:length(techs)
    % calcualate uniform thickness
    cas{j}         = subplot(1, techs_num, j);
    ca             = cas{j};
    imagesc(ca, xs{j}, zs{j}, Bscans_bg{j});
    % figure setup
    colormap(gray);
    xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
    % draw the line indicating the front surface
    hold on;
    plot(ca, xs{j}, front_line{j}, 'r-', 'linewidth', 2);
    % draw the line indicating the rear surface
    [row_inph, col_inph] = find(Bscans_track{j});
    row_inph             = row_inph  / fxs{j} * 1e3;
    col_inph             = col_inph / fss{j} * 1e6;
    hold on;
    scatter(ca, row_inph, col_inph, 8, 'cyano', 'd', 'filled');
    hold on;
    plot(ca, xs{j}, rear_line{j}, 'magenta-', 'linewidth', 2);
    hold on;
    set(ca, 'fontsize', 16, 'Fontname', 'times new Roman');
    set(ca, 'linewidth', 2);
    hold on;
    p = 1;
    while p <= 24
        scatter(xs{j}, front_line{j} + p * usTOF, 3, [0.8,0.8,0.8], 'o', 'filled');
        p = p + 1;
    end
    if j==1 % first col.
        legend({'Front-surface interply', 'Interplies excluding surfaces', 'Back-surface interply', 'Uniformly-spaced interplies'});
        ylabel(ca, 'TOF (\mus)', 'fontsize', 16, 'fontname', 'times new roman');
%         cl   = caxis;
        % adjust the pos.
        p    = get(ca, 'Position');
        p(1)    = p(1) - 0.07;
        p(2)    = p(2) + 0.02;
        p(3)    = p(3) + 0.08;
        set(ca, 'Position', p);
        h = colorbar;
        set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. amp. (arb.)');
    elseif j==3 % last col.
        p       = get(ca, 'Position');
        p(1)    = p(1) - 0.01;
        p(2)    = p(2) + 0.02;
        p(3)    = p(3) + 0.08;
        set(ca, 'Position', p);
%         caxis(cl);
        h = colorbar;
        set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. amp. (arb.)');
        set(ca, 'YTickLabel', []);
    else
%         caxis(cl);
        set(ca, 'YTickLabel', []);
        % adjust the pos.
        p       = get(ca, 'Position');
        p(1)    = p(1) - 0.04;
        p(2)    = p(2) + 0.02;
        p(3)    = p(3) + 0.08;
%         p(4)    = p(4) + 0.01;
        set(ca, 'Position', p);
        h = colorbar;
        set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. amp. (arb.)');
        set(ca, 'YTickLabel', []);
    end
    title(techs{j}, 'fontsize', 16, 'fontname', 'times new roman');
    xlabel({'\fontname{times new roman} X displacement (mm)', figlabel{j}}, 'fontsize', 16);
end

            
end

