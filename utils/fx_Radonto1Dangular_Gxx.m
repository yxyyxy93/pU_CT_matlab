function [R, anguler_1D, max_angle_I, xp, theta] = fx_Radonto1Dangular_Gxx(inam_C_scan, centers, radii)
% radon transform and reduce the 2D RDT iamge to 1D angular distribution
% select the max point in 1D angular distribution as well

% creat a round mask
mask = fx_createCirclesMask(inam_C_scan, centers, radii);

% use the mask to select a round region
inam_C_scan_mask = inam_C_scan .* mask;
x1 = centers(2) - radii;
x2 = centers(2) + radii;
y1 = centers(1) - radii;
y2 = centers(1) + radii;
inam_C_scan_mask = inam_C_scan_mask(x1: x2,y1: y2);

theta = 0:179;

% timer
tic;

[R, xp] = radon(inam_C_scan_mask, theta);

timeElapsed = toc;
disp(['radon: ', num2str(timeElapsed)]);

% calculating the gradient by central-different kernel
[Gx, ~] = imgradientxy(R, 'central');
[Gxx, ~] = imgradientxy(Gx, 'central');
anguler_1D = sum(abs(Gxx), 1);

timeElapsed = toc;
disp(['1D anguler: ', num2str(timeElapsed)]);

% calculating the maximum, for fft2
% anguler_1D = max(abs(R), [], 1);

[~, max_angle_I] = max(anguler_1D);

end

