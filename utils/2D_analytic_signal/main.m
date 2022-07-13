clc;
close all;
fclose all;
clear;


%% read image

% % filepath = 'C:\Users\xiayang\OneDrive - UGent\lab\cpp_work\2D_AS_v1\images';
% filepath = 'C:\Users\xiayang\OneDrive - UGent\lab\Optical_Microscopy_results\WovenFiber_optical\woven_patterns';
% 
% filename = strcat(filepath, '\woven_pattern2.png');
% % \original_i2d.png
% % \woven_pattern1
% % \woven_pattern2
% 
% img_origin = imread(filename);
% if size(img_origin,3)==3
%     img_gray = double(rgb2gray(img_origin));
% end
% 
% img_gray = max(img_gray) - img_gray;
%  - mean(img_gray, 'all')

%% construct i2D signal
[X,Y] = meshgrid(-17:0.05:17);

a       = 1;
theta_1 = 0;
theta_2 = 90;
phase  = 0;

img_gray = a*cos(X.*cosd(theta_1) + Y.*sind(theta_1) + phase); % i1D

img_gray = img_gray + a*cos(X.*cosd(theta_2) + Y.*sind(theta_2) + phase); %i1D + i1D = i2D

%%
clc;

s_0    = 10;
lambda = 0.5;
k      = 1;

s_c = s_0*lambda^(k-1); % coarse scale space parameter
s_f = s_c * lambda; % fine scale space parameter

fprintf("s_c: %0.2f \n", s_c);
fprintf("s_f: %0.2f \n", s_f);

% s_c = 10;  % coarse scale space parameter
% s_f = 4; % fine scale space parameter

% output relative bandwidth and center frequency
fprintf("relative bandwidth: %0.2f \n", lambda);
fprintf("octave number: %0.2f \n", k);
fprintf("wavelength: %0.2f \n", 2*pi*s_c*(lambda-1) / log(lambda));



%% 
mask_size = 5 * s_0;
[X_mesh, Y_mesh] = meshgrid(-mask_size: mask_size, -mask_size: mask_size);

kernel_1_f = fx_2DHilbertKernel_1st(X_mesh, Y_mesh, s_f);
kernel_1_c = fx_2DHilbertKernel_1st(X_mesh, Y_mesh, s_c);

kernel_2_f = fx_2DHilbertKernel_2nd(X_mesh, Y_mesh, s_f);
kernel_2_c = fx_2DHilbertKernel_2nd(X_mesh, Y_mesh, s_c);

fp_kernel   = kernel_1_f * s_f - kernel_1_c * s_c;
q1_x_kernel = (kernel_1_f - kernel_1_c) .* X_mesh;
q1_y_kernel = (kernel_1_f - kernel_1_c) .* Y_mesh;

q2_xx_kernel = (kernel_2_f - kernel_2_c) .* X_mesh.^2;
q2_xy_kernel = (kernel_2_f - kernel_2_c) .* X_mesh.*Y_mesh;
q2_yy_kernel = (kernel_2_f - kernel_2_c) .* Y_mesh.^2;

% fp_kernel   = kernel_1_f * s_f ;
% q1_x_kernel = (kernel_1_f) .* X_mesh;
% q1_y_kernel = (kernel_1_f) .* Y_mesh;
% 
% q2_xx_kernel = (kernel_2_f) .* X_mesh.^2;
% q2_xy_kernel = (kernel_2_f) .* X_mesh.*Y_mesh;
% q2_yy_kernel = (kernel_2_f) .* Y_mesh.^2;


% 
% 
% figure, subplot(1,3,1);
% surf(fp_kernel);
% colormap(gray);
% subplot(1,3,2);
% surf(q1_x_kernel);
% colormap(gray);
% subplot(1,3,3);
% surf(q1_y_kernel);
% colormap(gray);
% figure, subplot(1,3,1);
% surf(q2_xx_kernel);
% colormap(gray);
% subplot(1,3,2);
% surf(q2_xy_kernel);
% colormap(gray);
% subplot(1,3,3);
% surf(q2_yy_kernel);
% colormap(gray);
%
% 
f_s          = 1/2 * (q2_xx_kernel + q2_yy_kernel);
f_plus       = q2_xy_kernel;
f_plus_minus = 1/2 * (q2_xx_kernel - q2_yy_kernel);

figure, subplot(1,3,1);
surf(f_s);
colormap(gray);
subplot(1,3,2);
surf(f_plus);
colormap(gray);
subplot(1,3,3);
surf(f_plus_minus);
colormap(gray);


%% convolution 
[m, n]  = size(fp_kernel);
fftSize = size(img_gray) + [m, n] - 1;

f_p  = ifftn(fftn(img_gray, fftSize) .* fftn(fp_kernel, fftSize), 'symmetric');
f_x  = ifftn(fftn(img_gray, fftSize) .* fftn(q1_x_kernel, fftSize), 'symmetric');
f_y  = ifftn(fftn(img_gray, fftSize) .* fftn(q1_y_kernel, fftSize), 'symmetric');
f_xx = ifftn(fftn(img_gray, fftSize) .* fftn(q2_xx_kernel, fftSize), 'symmetric');
f_xy = ifftn(fftn(img_gray, fftSize) .* fftn(q2_xy_kernel, fftSize), 'symmetric');
f_yy = ifftn(fftn(img_gray, fftSize) .* fftn(q2_yy_kernel, fftSize), 'symmetric');

f_pm      = 0.5 * (f_xx - f_yy);
f_s       = 0.5 * f_p;
cos_alpha = sqrt(f_pm.^2 + f_xy.^2) ./ abs(f_s);
f_h       = sqrt((1+cos_alpha) ./ 2);
q         = (f_x.^2 + f_y.^2) .*2 ./ (1+cos_alpha);
%
phase = atan2(sqrt(q), f_p);
ampli = 0.5 * sqrt(f_p.^2 + q);
apexa = atan2(real(sqrt(f_s.^2 - f_xy.^2 - f_pm.^2)), sqrt(f_xy.^2 + f_pm.^2));
orien = 0.5 * atan2(f_xy, f_pm); % mean orientation
% orien    = atan2(f_y, f_x);
% mask_ori = phase<=0.1;
% orien(mask_ori) = 0.5 * atan2(f_xy(mask_ori), f_pm(mask_ori)) + pi/2;

% cut to shape
ex = 0;
f_p = f_p(m+1+ex:end-m-ex,     n+1+ex:end-n-ex);
phase = phase(m+1+ex:end-m-ex, n+1+ex:end-n-ex);
orien = orien(m+1+ex:end-m-ex, n+1+ex:end-n-ex);
ampli = ampli(m+1+ex:end-m-ex, n+1+ex:end-n-ex);
apexa = apexa(m+1+ex:end-m-ex, n+1+ex:end-n-ex);

% mask
mask = apexa <=0;

% % surf plot
% figure,
% subplot(2,3,1), imagesc(img_gray); colorbar;
% title('original');
% subplot(2,3,2);
% h =surf(f_p); colorbar;
% set(h,'LineStyle','none');
% title('signal in possion scale space');
% subplot(2,3,3);
% h = surf(phase); colorbar;
% set(h,'LineStyle','none')
% title('phase');
% orien_image = orien;
% % orien_image(~mask) = nan;
% ca = subplot(2,3,4);
% h = surf(orien_image); colorbar;
% set(h,'AlphaData',~isnan(orien_image));
% set(h,'LineStyle','none');
% colormap(ca, hsv); colorbar;
% title('orientation');
% subplot(2,3,5);
% h = surf(ampli); colorbar;
% set(h,'LineStyle','none');
% title('amplitude');
% subplot(2,3,6);
% h = surf(apexa); 
% set(h,'LineStyle','none');
% colorbar;
% title('apex angle');

% imagesc plot
figure,
subplot(2,3,1), h = imagesc(img_gray); colorbar;
set(h,'AlphaData',~isnan(img_gray));

title('original');
subplot(2,3,2);
h =surf(f_p);
set(h,'LineStyle','none');
colorbar;
title('signal in possion scale space');
subplot(2,3,3);
h = imagesc(phase); colorbar;
title('phase');
orien_image = orien;
orien_image(~mask) = nan;
ca = subplot(2,3,4);
h = imagesc(orien_image); colorbar;
set(h,'AlphaData',~isnan(orien_image));
colormap(ca, hsv); colorbar;
title('orientation');
subplot(2,3,5);
h = imagesc(ampli); colorbar;
title('amplitude');
subplot(2,3,6);
h = imagesc(apexa); 
colorbar;
title('apex angle');

