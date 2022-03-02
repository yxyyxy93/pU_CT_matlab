function cf = fx_errorsbars_tightplot(t_space, snrs, techs, errors, stds, std_mean_threshold)
% plot 9 errors bar in 3*3 figure tightly
% snrs:     3 levels of snr
% techs:   3 kinds of techs
% errors: 3*3 cell of errors
% stds:    3*3 cell of stds

% plot in errorbars
cf        = figure('Name', ['best_results_' 'errorbars']);
set(cf, 'Position', [0, 0, 1600, 800], 'color', 'white');
cas       = cell(length(techs)*length(snrs));
snr_num   = length(snrs);
techs_num = length(techs);   
alpha     = 0.3;
for i = 1:snr_num
    for j = 1:length(techs)
        cas{i*snr_num+j-snr_num} = subplot(snr_num, techs_num, i*snr_num+j-snr_num);
        ca                       = cas{i*snr_num+j-snr_num};
        %
        Origin_mean = errors{i, j};
        Origin_std  = stds{i, j};
        untrackable = Origin_mean>std_mean_threshold(1) | Origin_std>std_mean_threshold(2);
        h           = errorbar(ca, t_space(untrackable), Origin_mean(untrackable)*100, ...
            Origin_std(untrackable)*100, 's', 'MarkerSize',5, 'MarkerFaceColor', [1 1 1]*(1-alpha), ...
            'CapSize', 10, 'LineWidth', 2, 'LineStyle', 'none');
        h.Color     = [1 1 1]*(1-alpha);
        % Set transparency
        set([h.Bar, h.Line], 'ColorType', 'truecoloralpha', ...
            'ColorData', [h.Line.ColorData(1:3); 255*alpha]);
        hold on;
        h           = errorbar(ca, t_space(~untrackable), Origin_mean(~untrackable)*100, ...
            Origin_std(~untrackable)*100, 's', 'MarkerSize', 5, 'MarkerEdgeColor', 'red', ...
            'MarkerFaceColor','red', 'CapSize', 10, 'LineWidth', 2, 'LineStyle', 'none');
        h.Color     = [0.00,0.45,0.74];
%         %
%         errorbar(cas{i*snr_num+j-snr_num}, t_space, errors{i, j}*100, stds{i, j}*100, '-s',...
%             'MarkerSize', 5, 'Linewidth', 2,'MarkerEdgeColor','red','MarkerFaceColor','red');
%         ytickformat('percentage');
        set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
        set(ca, 'linewidth', 2);
        if j==1
            ylabel({['SNR=' num2str(snrs(i))], 'Errors (\times 100%)'}, 'fontsize', 16, 'fontname', 'times new roman');
        else
            set(ca, 'YTickLabel', []);
        end
        if i==1
            title(techs{j}, 'fontsize', 16, 'fontname', 'times new roman');
            set(ca, 'XTickLabel', []);
        elseif i==3
            xlabel({'\fontname{times new roman} TOF (\mus)'}, 'fontsize', 16);
        else
            set(ca, 'XTickLabel', []);
        end
        grid on;
        p    = get(ca, 'Position');
        % tight the subplot
        p(3) = p(3) + 0.06;
        p(4) = p(4) + 0.05;
        set(ca, 'Position', p);
        % set ylim;
        yl = ylim;
        ylim(yl/3);
    end
end

for i = 1:length(cas)-1
    linkaxes([cas{i}, cas{i+1}], 'xy');
end