function fig = fx_showBS_one(fignum, inph, index, front_I, rear_I, B_type, name_colormap, titlestr, x, y)
% show b scan

switch nargin
    case 7
        titlestr = 'B scan section of ';
        B_typestr = B_type;
    case 8
        disp(titlestr);
        B_typestr = '';
end

% save and output the handles
fig  = figure(fignum);

if (B_type == 'x')
    % B scan image
    B_scan = squeeze(inph(:, index, :));
    % front line
    front_line = squeeze(front_I(:, index));
    % rear line
    rear_line = squeeze(rear_I(:, index));
elseif (B_type == 'y')
    % B scan image
    B_scan = squeeze(inph(index, :, :));
    % front line
    front_line = squeeze(front_I(index, :));
    % rear line
    rear_line = squeeze(rear_I( index, :));
end
imagesc(y, x, B_scan');
hold on;
% figure setup
title([B_type ' = ', num2str(index)], 'fontsize', 16, 'Fontname', 'times new Roman');
colormap(name_colormap);
h = colorbar;
set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16}');
set(gca, 'Fontname', 'times new Roman');
xlabel('\fontname {times new roman}', 'fontsize', 16);
ylabel('\fontname {times new roman}', 'fontsize', 16);
% draw the line indicating the front surface
hold on
plot(front_line, 'red--', 'linewidth', 3);
% draw the line indicating the rear surface
hold on
plot(rear_line, 'red--', 'linewidth', 3);

set(gcf, 'Position', [0, 0, 600, 400], 'color', 'white');
set(gca, 'fontsize', 16);
set(gca, 'fontname', 'Times New Roman');
set(gca, 'linewidth', 1.5);

end

