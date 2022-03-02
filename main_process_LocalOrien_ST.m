% read the file
clc;
close all;
fclose all;
clear;

%% Define the class
% read the preprocessed .mat
[Filename1,Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename              = strcat(Pathname1, Filename1);
process               = class_process_RFdata(filename);

% read the settings
[Filename1,Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
settinig_file         = strcat(Pathname1, Filename1);
process               = process.loadSettings(settinig_file);

% read the data from the struct, or reset the dataset
process               = process.read_origin_data; % read (reset) the dataset
x                     = 10;
y                     = 10;

%% EDA
process.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500];
% window  = [21 280 121 380 1 1500]; 
window  = [101 300 201 400 1 1500];  % no damage
x_step  = 2;
y_step  = 2;
process = process.cut_edges(window, x_step, y_step);
process.show_img_shape;

% check infq
s       = squeeze(process.img(x, y, :));
[infq, inph, inam, hs]= fx_ht_inFqPhAm(s', process.fs);
figure('Name', 'Inst. freq.');
plot(infq);

process.show_hilbert_Ascan(x, y);

xslice = 50 / process.fx * 1e3;
yslice = 50 / process.fy * 1e3;
zslice = 1000 / process.fs * 1e6;
process.show_inph_3d('img_hil', xslice, yslice, zslice);
process.show_unwraped_inph_3d('img_hil', xslice, yslice, zslice);

% filter and demo Scaleogram
x     = 20;
y     = 20;
f0    = linspace(1e6, 15e6, 100);
sigma = 0.7;
process.show_logGabor_Scaleogram(f0, sigma, x, y);

% 3d display
fx_Scrollable_3d_view(abs(process.img_hil_filter));

%% read reference signal
[Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
filename               = strcat(Pathname1, Filename1);
x                      = 10;
y                      = 10;
process                = process.read_refer(filename, x, y);
% cut the ref. signal
win                    = [0.1e-6 2e-6]; % unit: s
process                = process.cut_reference_signal(win);

% deconvolution
q_factor    = 1e-2;
ker_wid     = 200;
process = process.apply_deconvolutions_onlyWiener...
    (q_factor, ker_wid);

%% surface calculation
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
MinPD        = 25;
MinPH        = 0.2;  % these 2 parameters need to be changed for surface estimation.
PropertyName = 'img_hil';
process.show_Ascan_inam_peaks(77, 77, MinPD, MinPH, PropertyName); % x, y

MinPH        = 1.5e-3;  % these 2 parameters need to be changed for surface estimation.
PropertyName = 'img_WienerDeconv';
process.show_Ascan_inam_peaks(19, 78, MinPD, MinPH, PropertyName); % x, y

% surface searching
max_len = 1000;
process = process.find_front_amp(MinPD, MinPH, PropertyName, max_len);
%
filename_fig = filename;
process.show_surfaces(filename_fig(1:end-4));
%
process      = process.smooth_rear_I;
process.show_surfaces(filename_fig(1:end-4));
%
process      = process.recover_surface;

%% 3d ply track
% logGabor filtered
f0      = 6.3e6;
sigma   = 0.75;
% Ascan demo
process.demo_logGabor_plytrack(x, y, f0, sigma);
% 
process = process.Filter_logGabor(f0, sigma, 'img');
process = process.track_interply('img_hil_filter');
xslice  = 50 / process.fx * 1e3;
yslice  = 100 / process.fy * 1e3;
zslice  = [];

% process.show_track_interply(xslice, yslice, zslice);
% process.show_B_scan_interply('x', 80, 1:200, 'img_hil_filter');

% 3d interply track
process = process.show_interply_track_3D_knn(xslice, yslice, zslice);

%%
% extract the scatter of the interply track with the markers indicating the layers.
% convert to 3d scattering points and initial the props for clustering
process = process.show_interply_track_3D;

% KNN search
K       = 400;
P       = 2; % Euclidean distance
nbins   = 100;
process = process.knn_search(K, P, nbins);

% Clustering by KNN 
threshold_dis = 22;
process = process.knn_clustering(threshold_dis);

% release the memory
process = process.clear_knn_properties;

% show one layer interply track in C-scan mode
process.show_one_interply_knn(23, 3000);

% 3d interply track
process = process.show_interply_track_3D_knn(xslice, yslice, zslice);

%% structrul tensor 
% out-of-plane angles
% logGabor filtered
f0      = 7.5e6;
sigma   = 0.2;
% Ascan demo
process.demo_logGabor_plytrack(x, y, f0, sigma);
process = process.Filter_logGabor(f0, sigma, 'img');

% smoothing scale
sigma1       =[3, 3, 3];
% integration scale
sigma2       = [3, 3, 3];
PropertyName = 'img_hil_filter';
process      = process.structural_tensor(sigma1, sigma2, PropertyName);

% display the angles
medf_kernel = [3, 3, 3];
xslice      = 50 / process.fx * 1e3;
yslice      = 50 / process.fy * 1e3;
zslice      = [];
process.show_angles_ST(medf_kernel, xslice, yslice, zslice);

%% local in-plane angles extract one image demo 
% 2D demo
% total_plies = 24;
% [angle, process] = process.show_one_Radontransform_vitualprofile(5, 0.5, 'img_hil_filter', center, radii, total_plies);
% Radontransform
% define the C_scan_inam as well
% PropertyName = 'img_hil';
PropertyName = 'img';
ply          = 22;
ratio        = 0.5;
process      = process.define_plywise_inamCscan(ply, ratio, PropertyName);
% curvelet denoising
sigma        = 5e-4;
imagename    = 'C_scan_inam_plywise';
process      = process.curvelet_denoise_inamCscan(sigma, imagename);

imagename = 'C_scan_inam_denoise';
theta   = 0:1:180;
ply_wise_Cscan_inam = cell(1, 24);
anguler_1D_ply_wise = cell(1, 24);
for ply = 1:24
    ratios  = 0.2:0.1:0.8;
    process = process.define_plywise_inamCscan(ply, ratios, 'img_hil_filter');
    % curvelet denoising
    sigma   = 1e-2;
    process = process.curvelet_denoise_inamCscan(sigma, imagename);
    [R, anguler_1D, max_angle_I, xp] = fx_Radonto1Dangular...
        (process.C_scan_inam_denoise, center, radii, theta);
    % plot the 1D anguler (sum of gradient by central-different kernel)
    figure, plot(theta, anguler_1D, 'linewidth', 2);
    anguler_1D_ply_wise{ply} = anguler_1D;
    ply_wise_Cscan_inam{ply} = process.C_scan_inam_denoise;
end

% display all images in on figure
close all;
figure('Name', 'ply_wise_Cscan_inam')
ply_wise_Cscan_inam_bg = 0; % the background
for ply = 1:24
    subplot(6, 4, ply)
    imagesc(ply_wise_Cscan_inam{ply}); colormap gray; title(num2str(ply));
    ply_wise_Cscan_inam_bg = ply_wise_Cscan_inam_bg + ply_wise_Cscan_inam{ply};
end
%
figure('Name', 'ply_wise_anguler_1D')
for ply = 1:24
    subplot(6, 4, ply)
    plot(theta, anguler_1D_ply_wise{ply}, 'linewidth', 2);
end

% PropertyName = 'img_hil_filter';
% dis_to_front = 3.7;
% process      = process.define_parallel_inamCscan...
%     (dis_to_front, PropertyName, center, radii);

% curvelet denoising
sigma   = 1e-2;
process = process.curvelet_denoise_inamCscan(sigma);

% RT
theta     = 0:5:180;
radiis    = 20;
imagename = 'C_scan_inam_denoise';
process   = process.compute_orientation_by_RT_correct(radiis, theta, imagename);
process.show_orientation_by_ID_RT(radiis, theta);
imagename = 'C_scan_inam';
process   = process.compute_orientation_by_RT_correct(radiis, theta, imagename);
process.show_orientation_by_ID_RT(radiis, theta);

% 2d log-Gabor fitler
wavelength  = [2 4 8 16];
orientation = 0:5:180;
% process = process.show_logGabor_filter('img_hil_filter', wavelength, orientation);
SFB         = 1; % [0.5 2.5]
SAR         = 0.5; % [0.23 0.92]

% K: controls how much smoothing is applied to the Gabor magnitude responses.
K = 1e-2;
imagename = 'C_scan_inam_denoise';
process     = process.compute_logGabor_filter_withoutFig('img_hil_filter', wavelength, orientation, SFB, SAR, imagename);
process.show_orientation_by_ID_allwl(wavelength, orientation, K);
imagename = 'C_scan_inam';
process     = process.compute_logGabor_filter_withoutFig('img_hil_filter', wavelength, orientation, SFB, SAR, imagename);
process.show_orientation_by_ID_allwl(wavelength, orientation, K);

%% local in-plane angles extract
% 3D volumetric 
% logGabor filtered
f0      = 7.5e6;
sigma   = 0.5;
% Ascan demo
process.demo_logGabor_plytrack(x, y, f0, sigma);
% 
process = process.Filter_logGabor(f0, sigma, 'img');

% inplane orientation extraction for 3D image
ratios      = 0.1:0.1:0.9;
total_plies = 24;
SFB         = 1; % [0.5 2.5]
SAR         = 0.5; % [0.23 0.92]
K           = 0.01;
wavelength  = [4 8 16 32];
orientation = 0:5:180;
process     = process.extract_local_orientation_3D('img_hil_filter', wavelength, orientation, ratios, total_plies, SFB, SAR, K);

% show slice
xslice      = 50 / process.fx * 1e3;
yslice      = 50 / process.fy * 1e3;
zslice      = [];
process.show_inplane_direction_3D(xslice, yslice, zslice);

% % inplane orientation extraction by LG (dis_to_front)for 3D image
% distances_to_front = 0.1:0.1:3.7;
% process            = process.extract_local_orientation_3D_parallel_allwl...
%     ('img_hil_filter', wavelength, orientation, distances_to_front);
% % show slice
% xslice             = 100 / process.fx * 1e3;
% yslice             = 100 / process.fy * 1e3;
% zslice             = [];
process.show_inplane_direction_3D_ID(xslice, yslice, zslice);

[Filename1,Pathname1] = uigetfile({'*.mat'}, 'select the file');
loadname              = strcat(Pathname1, Filename1);
process               = process.load_local_orientaion_3D(loadname);
xslice                = 10;
yslice                = 10;
process.show_inplane_direction_3D_ID(xslice, yslice, zslice);

%% distance to front

% RT
% PropertyName       = 'img_WienerDeconv'; % deconvolved signal
PropertyName       = 'img_hil_filter'; % orginal signal
radius             = 20;
theta              = 0:5:180; % degree
distances_to_front = 0.05:0.05:3.75;
process            = process.extract_local_orientation_RT_3D_parallel(...
    PropertyName, radius, theta, distances_to_front);

% 2D log-Gabor filter
wavelength         = [2 4 8 16 32];
orientation        = 0:1:180;
distances_to_front = 0.01:0.01:3.75;
process            = process.extract_local_orientation_3D_parallel_allwl(...
    PropertyName, wavelength, orientation, distances_to_front);

% show slice
xslice      = 50 / process.fx * 1e3;
yslice      = 50 / process.fy * 1e3;
zslice      = [];
process.show_inplane_direction_3D_ID(xslice, yslice, zslice);
            