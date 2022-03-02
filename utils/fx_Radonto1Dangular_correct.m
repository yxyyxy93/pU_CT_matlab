function [R, anguler_1D, max_angle_I, xp] = fx_Radonto1Dangular_correct(inam_C_scan, centers, radii, theta)
% radon transform and reduce the 2D RDT iamge to 1D angular distribution
% select the max point in 1D angular distribution as well
% theta: the theta space

x1 = centers(2) - radii;
x2 = centers(2) + radii;
y1 = centers(1) - radii;
y2 = centers(1) + radii;
% % creat a round mask
% mask = fx_createCirclesMask(inam_C_scan, centers, radii);
% % use the mask to select a round region
% inam_C_scan_mask = double(inam_C_scan) .* mask;
% inam_C_scan_mask = inam_C_scan_mask(x1: x2,y1: y2);

inam_C_scan_mask = inam_C_scan(x1: x2, y1: y2);
[xx, yy]         = meshgrid(-radii:radii, -radii:radii);
mask             = (hypot(xx, yy)<=radii);
inam_C_scan_mask = double(inam_C_scan_mask) .* mask;

% timer
[R, xp] = radon(inam_C_scan_mask, theta);

% cut the edge of the radon image
R_win  = R(ceil(end / 2 - radii + 1): floor(end / 2 + radii), :);
xp_win = xp(ceil(end / 2 - radii + 1): floor(end / 2 + radii));

% consider the change of length of the intergration line 
integration_line_length =  1 .* 2 .* (radii.^2 - xp_win.^2).^0.5;
integration_line_length = repmat(integration_line_length, 1, size(R_win, 2));
R_divide                = R_win ./ integration_line_length;
% % without considering the change of length of the intergration line 
% R_divide                = R_win;

% Gyy = diff(Gy, 1, 1); % 2nd
% takes the sum of the absolute value of the first derivative of the Radon transform in the radial
% Gy         = diff(R_divide, 1, 1);
%!this line produce strange results! Gy = gradient(R_divide, 1); % using a central-difference kernel 
Gy         = (R_divide(3:end, :) - R_divide(1:end-2, :)) / 2; % using a central-difference kernel 

Gy_abs     = abs(Gy);
anguler_1D = sum(Gy_abs, 1);

% anguler_1D = std(R_divide, 0, 1) ./ mean(R_divide, 1);
[~, max_angle_I] = max(anguler_1D);


end

