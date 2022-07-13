clc;
close all;
fclose all;
clear;

%% **************************************
% Define the class
% read the preprocessed .mat
[Filename1,Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename              = strcat(Pathname1, Filename1);
process3              = class_process_RFdata(filename);

% read the settings
[Filename1,Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
settinig_file         = strcat(Pathname1, Filename1);
process3              = process3.loadSettings(settinig_file);

% read the data from the struct, or reset the dataset
process3              = process3.read_origin_data; % read (reset) the dataset
% 

%% ********** surface (8p waviness sample)
% 2D demo
% total_plies = 24;
% [angle, process] = process.show_one_Radontransform_vitualprofile(5, 0.5, 'img_hil_filter', center, radii, total_plies);

% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
MinPD        = 50;
MinPH        = 0.13;  % these 2 parameters need to be changed for surface estimation.
% PropertyName = 'img_WienerDeconv'; bh
PropertyName = 'img_hil';
process3.show_Ascan_inam_peaks(210, 304, MinPD, MinPH, PropertyName); % x, y

% surface calculation
max_len                = 700;
threshold_delamination = 1;
process3               = process3.find_front_amp(...
    MinPD, MinPH, PropertyName, max_len, threshold_delamination);
filename_fig           = filename;
process3.show_surfaces(filename_fig(1:end-5));

%
process3.smooth_rear_I;
process3.show_surfaces(filename_fig(1:end-5));

%% window cutting
x_step   = 1;
y_step   = 1;
window   = [1 350 1 550 1 1000];
% x_step                = 1;
% y_step                = 1;
process3 = process3.cut_edges(window, x_step, y_step);

process3.show_img_shape;

%% ********** one-plane EDA
% define the C_scan_inam as well
% PropertyName = 'img_WienerDeconv';
f0       = 15e6;
sigma    = 0.7; 
process3 = process3.Filter_logGabor(f0, sigma, 'img');

PropertyName = 'img_hil_filter';
% PropertyName = 'img_hil';
ratio        = 3.1/8;
process3     = process3.define_parallel_inamCscan...
    (ratio, PropertyName);
%
imagename = 'C_scan_inam';
% curvelet denoising
sigma     = 5e-4;
process3  = process3.curvelet_denoise_inamCscan(sigma, imagename);

% Radontransform
theta         = 1:1:180;
radiis        = 10;
angle_compens = 0;
process3      = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);

% 2d log-Gabor fitler
wavelength  = 4:2:20;
orientation = 1:1:180;
SFB         = 1; % [0.5 2.5]
SAR         = 0.7; % [0.23 0.92]

% imagename = 'C_scan_inam_denoise';
% process3  = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% process3.show_orientation_by_ID_allwl(wavelength, orientation, K);
imagename     = 'C_scan_inam';
process3      = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K             = 1e0;
angle_compens = 0;
process3.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);

%% 3D information Diagram
PropertyName  = 'img_hil';

wavelength  = 8:1:12;
orientation = 0:1:180; % degree
SFB         = 2;       % [0.5 2.5]
SAR         = 0.5;     % [0.23 0.92]
z_theta     = 0:5:60;
K           = 1;   
z_range     = [1 1200];
%
tic;
process3 = process3.extract_local_orientation_3DID...
    (PropertyName, wavelength, orientation, z_theta, z_range, K);
toc;

% % 3D ID plane
% y_ori = 0:1:20;
% z_ori = 0:1:20;
% 
% tic;
% process3 = process3.extract_local_orientation_3DID_plane...
%     (PropertyName, wavelength, y_ori, z_ori, z_range, K);
% toc;

% show slice
xslice        = 50 / process3.fx * 1e3;
yslice        = 50 / process3.fy * 1e3;
zslice        = [];
mfsize        = [1 1 1];
angle_compens = 14;
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

process3.statistic_angular_distribution(angle_compens);


%% one-plane EDA 2D analytic-signal
% PropertyName = 'img_hil_filter';
PropertyName = 'img_hil';
% PropertyName = 'img_WienerDeconv';
% PropertyName = 'img_WienerDeconv_AR'';

process3 = process3.define_parallel_inamCscan(2.5/8, PropertyName);
imagename   = 'C_scan_inam';
clc;

s_0    = 1;
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

mask_size = 10;
process3.compute_2D_analytic_signal_image(PropertyName, s_f, s_c, mask_size, imagename);

% % scale adapted
% scale_map = ones(size(process3.(imagename))- 2*mask_size-2);
% ampli_map = zeros(size(process3.(imagename)) - 2*mask_size-2);
% orien_map = zeros(size(process3.(imagename)) - 2*mask_size-2);
% apexa_map = zeros(size(process3.(imagename)) - 2*mask_size-2);
% for s_0 = 0.5:10
%     s_c = s_0*lambda^(k-1); % coarse scale space parameter
%     s_f = s_c * lambda; % fine scale space parameter
%     [~, orien, ampli, apexa] = process3.compute_2D_analytic_signal_image(PropertyName, s_f, s_c, mask_size, imagename);
%     mask = ampli > ampli_map;
%     ampli_map = mask .* ampli + (1-mask) .* ampli_map;
%     orien_map = mask .* orien + (1-mask) .* orien_map;
%     apexa_map = mask .* apexa + (1-mask) .* apexa_map;
%     scale_map(mask) = s_0;
% end
% close all;
% orien_map(apexa_map>=1) = nan;

%% 'distance to front and rear' in-plane orientation extraction
% RT
% PropertyName     = 'img_WienerDeconv'; % deconvolved signal
% PropertyName  = 'img_hil_filter'; % orginal signal
PropertyName  = 'img_hil';

radius        = 10;
theta         = 1:1:180; % degree
ratio         = 0/24:0.1/24:24/24;
sigma_denoise = 0;
tic;
process3      = process3.extract_local_orientation_RT_3D_parallel(...
    PropertyName, radius, theta, ratio, sigma_denoise);
toc;

% show slice
xslice        = 150 / process3.fx * 1e3;
yslice        = 50 / process3.fy * 1e3;
zslice        = []; % us
angle_compens = -7;
process3.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);

% 2D log-Gabor filter
wavelength         = 4:2:20;
orientation        = 1:1:180;
ratio              = 0/24:1/24:24/24;
K                  = 0.5e0;
% sigma              = 5e-4;
sigma              = 0;
tic;
process3           = process3.extract_local_orientation_3D_parallel_allwl(...
    PropertyName, wavelength, orientation, ratio, K, sigma);
toc;

% show slice
xslice        = 50 / process3.fx * 1e3;
yslice        = 50 / process3.fy * 1e3;
zslice        = []; % us
angle_compens = 0;
mfsize        = [1 1 1];
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

%% *************** 3d ply track
% add noise to img.
% process = process.addnoise(25);
% origin without filter
threshold  = 0.01;  % max(inam) * threshold for distinguishing the back surface

% % process  = process.track_interply('img_hil');
% process3   = process3.track_interply_inph(threshold, 'img_hil');
% 
% %
% % process.show_track_interply(xslice, yslice, zslice);
% process3.show_oneinterply(5, 'nofilter', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply(15, 'nofilter', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply(23, 'nofilter', 3000, win_x, win_y, TOF_oneply);
% process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

% % filtered
f0        = 5.5e6;
sigma     = 0.7;

% for debug
% threshold = 0.05;
x         = 100;
y         = 98;
% process3.demo_logGabor_plytrack_inph_v2(x, y, f0, sigma, threshold);
process3.show_hilbert_Ascan(x, y);
process3.show_logGabor_Scaleogram(3e6:0.1e6:16e6, sigma, x, y);

% 0.676 us - 2.13 us

%  
process3 = process3.Filter_logGabor(f0, sigma, 'img');
% nol      = 17;
% process3 = process3.track_interply_2ndharmonic('img_hil_filter', nol);
process3 = process3.track_interply('img_hil_filter');

% use 2nd-harmonic for first and the last interplies, 
% use fundamental resonance for another interplies
f0_1     = 11e6;
sigma0_1 = 0.8;
f0_2     = 5.5e6;
sigma0_2 = 0.7;
nol      = 9;
process3 = process3.track_interply_hybrid(...
    'img', f0_1, sigma0_1, f0_2, sigma0_2, nol);

%%
B_type     = 'x';
index      = 200;
Bwin       = 1:200;
TOF_oneply = TOF_walls / 24; % 24 plies
process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

% 3d slices
xslice     = 175 / process3.fx * 1e3;
yslice     = 175 / process3.fy * 1e3;
zslice     = [];
process3.show_interply_track_3D_knn(xslice, yslice, zslice);
% process.show_track_interply(xslice, yslice, zslice);
win_x      = 1:80;
win_y       = 1:80;
process3.show_oneinterply(05, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(15, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(22, 'logGabor', 3000, win_x, win_y, TOF_oneply);

std_mean_threshold = [1, 0.22]; % modify here 
figname            = 'logGabor';
process3.thickness_estimation_v1(5.5, 24, win_x, win_y, figname);

% plot with uniformly spaced locaitons
threshold = 0.05;
nol       = 25;
TOF       = process3.track_show_B_scan_interply(...
    B_type, 80, Bwin, 'img_hil_filter', threshold, nol, std_mean_threshold);
% 
% % pos. err. calculation 
% % win_x      = 50:250;
% % win_y      = 50:250;
% layers_num = 25;
% process3   = process3.error_pos_calculate(win_x, win_y, layers_num);

% TOF    = process3.track_show_B_scan_interply(B_type, index, Bwin, 'img_hil_filter', threshold);
save TOF.mat TOF_walls;

%% one-plane EDA
% 2D demo
% define the C_scan_inam as well
% % filtered
f0       = 15e6;
sigma    = 0.7; 
process3 = process3.Filter_logGabor(f0, sigma, 'img');
%
PropertyName = 'img_hil_filter';
% PropertyName = 'img_hil';
% PropertyName = 'img_WienerDeconv';
ply          = 3;
ratio        = 0.5;
process3     = process3.define_plywise_inamCscan(ply, ratio, PropertyName);
% curvelet denoising
sigma        = 5e-3;
%
imagename     = 'C_scan_inam_plywise';
process3      = process3.curvelet_denoise_inamCscan(sigma, imagename);
% RT
theta         = 1:1:180;
radiis        = 10;
% imagename     = 'C_scan_inam_plywise';
% imagename     = 'C_scan_inam_denoise';
process3      = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
angle_compens = 0;
process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);

% 2d log-Gabor fitler
wavelength  = 4:1:20;
orientation = 1:1:180;
SFB         = 1; % [0.5 2.5]
SAR         = 0.7; % [0.23 0.92]

% K: controls how much smoothing is applied to the Gabor magnitude responses.
K             = 0.5e0;
process3      = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
angle_compens = 0;
process3.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);

%%
FolderName  = "F:\Xiayang\results\image_processing\";
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '12022021_11h25m10s.mat'));

% 10022021_15h49m02s       % waviness sample
% 12032021_11h14m41s   
%% save all figures 3rd paper
FolderName  = "F:\Xiayang\results\image_processing";   % the destination folder
FigList     = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  disp(fullfile(FolderName, FigName));
  set(0, 'CurrentFigure', FigHandle);
  saveas(gcf, strcat(FolderName, '\', FigName), 'epsc');
  saveas(gcf, strcat(FolderName, '\', FigName), 'pdf');
  saveas(gcf, strcat(FolderName, '\', FigName), 'fig');
end
