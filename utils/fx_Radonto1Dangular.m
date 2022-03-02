function [R, anguler_1D, max_angle_I, xp] = fx_Radonto1Dangular(inam_C_scan, centers, radii, theta)
% radon transform and reduce the 2D RDT iamge to 1D angular distribution
% select the max point in 1D angular distribution as well
% theta: the theta space
% creat a round mask
mask = fx_createCirclesMask(inam_C_scan, centers, radii);

% use the mask to select a round region
inam_C_scan_mask = double(inam_C_scan) .* mask;

x1               = centers(2) - radii;
x2               = centers(2) + radii;
y1               = centers(1) - radii;
y2               = centers(1) + radii;
inam_C_scan_mask = inam_C_scan_mask(x1: x2,y1: y2);

% timer
% tic;
[R, xp] = radon(inam_C_scan_mask, theta);
% % consider the change of length of the intergration line 
% integration_line_length = 2 * (radii.^2 - xp.^2).^0.5;
% integration_line_length = 1 * repmat(integration_line_length, 1, size(R, 2));
% R_divide = R ./ integration_line_length;

% cut the edge of the
R_win = R(ceil(end / 2 - radii + 1): floor(end / 2 + radii), :);
% % consider the change of length of the intergration line 
% xp_win = xp(ceil(end / 2 - radii + 1): floor(end / 2 + radii));
% integration_line_length =  1 .* 2 .* (radii.^2 - xp_win.^2).^0.5;
% integration_line_length = repmat(integration_line_length, 1, size(R_win, 2));
% R_divide = R_win ./ integration_line_length;
R_divide = R_win;

% timeElapsed = toc;
% disp(['radon: ', num2str(timeElapsed)]);

% calculating the gradient by central-different kernel
% [~, Gy] = imgradientxy(R_win, 'central');
Gy = diff(R_divide, 1, 1);

% Gyy = diff(Gy, 1, 1); % 2nd
Gy_abs = abs(Gy);
anguler_1D = mean(Gy_abs, 1) ./ mean(R_divide, 1);

% disp(num2str(mean(anguler_1D)));

% anguler_1D_gyy = mean(abs(Gyy), 1);
% [Gxy, Gyy] = imgradientxy(Gy, 'central');
% anguler_1D_Gyy = sum(abs(Gyy), 1);
[~, max_angle_I] = max(anguler_1D);


end

