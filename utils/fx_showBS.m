function [fig, ax]= fx_showBS(fignum, inph, num_bs,  step_bs, front_I, rear_I, B_type, name_colormap, titlestr)
% show b scan

switch nargin
    case 8
        titlestr = 'B scan section of '; 
        B_typestr = B_type;
    case 9
        disp(titlestr);    
        B_typestr = '';
end

% save and output the handles
fig = figure(fignum);
ax = cell(num_bs);

for i = 1: num_bs
    ax{i} = subplot(2, num_bs/2 , i);
    B_sec = i*step_bs;                  % the place of the section
    if (B_type == 'x')
        % B scan image
        B_scan = squeeze(inph(:, B_sec, :));
        % front line
        front_line = squeeze(front_I(:, B_sec));
        % rear line
        rear_line = squeeze(rear_I(:, B_sec));
    elseif (B_type == 'y')
        % B scan image
        B_scan = squeeze(inph(B_sec, :, :));
        % front line
        front_line = squeeze(front_I(B_sec, :));
        % rear line
        rear_line = squeeze(rear_I( B_sec, :));
    end
    imagesc(B_scan');
    hold on;
    % figure setup
    title([B_type ' = ', num2str(B_sec)], 'fontsize', 16, 'Fontname', 'times new Roman');   
    colormap(name_colormap);
    h = colorbar;
    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16}');
    set(gca, 'Fontname', 'times new Roman');
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    xlabel('\fontname {times new roman}', 'fontsize', 16);
    ylabel('\fontname {times new roman}', 'fontsize', 16);
    % draw the line indicating the front surface
    hold on
    plot(front_line, 'red.', 'markersize', 2);
    % draw the line indicating the rear surface
    hold on
    plot(rear_line, 'red.', 'markersize', 2);
end
set(gcf, 'Position', [0, 0, 800, 600], 'color', 'white');
suptitle([titlestr, B_typestr]);

end

