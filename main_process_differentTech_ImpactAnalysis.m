% read the file
clc;
close all;
fclose all;
clear;

%%
% 15 MHz Analytic-signal analysis
% Define the class
% read the preprocessed .mat
[Filename1, Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename               = strcat(Pathname1, Filename1);
process3               = class_process_RFdata(filename);

% read the settings
[Filename1, Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
settinig_file          = strcat(Pathname1, Filename1);

% read the data from the struct, or reset the dataset
process3               = process3.loadSettings(settinig_file);
process3               = process3.read_origin_data; % read (reset) the dataset

x                      = 10;
y                      = 10;

%% window cutting
window = [1 321 51 400 1 1500];
x_step = 1;
y_step = 1;
% window   = [11 101 11 201 100 1400];
x_step   = 1;
y_step   = 2;
window   = [1 350 51 750 100 1400];
% x_step                = 1;
% y_step                = 1;
process3 = process3.cut_edges(window, x_step, y_step);

process3.show_img_shape;

% shift the A scan in the time domain
process3 = process3.shift_A_scan(0);

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
MinPD        = 25;
MinPH        = 0.03;  % these 2 parameters need to be changed for surface estimation.
PropertyName = 'img_hil';
process3.show_Ascan_inam_peaks(78, 78, MinPD, MinPH, PropertyName); % x, y
process3.show_Ascan_inam_peaks(46, 69, MinPD, MinPH, PropertyName); % x, y

% surface calculation
max_len  = 1300;
% % the attenuation coefficient for defining the echo from delamiation.
% Ab         = 0.075;
% A0         = 0.83;
% h          = 1156 - 240; % data point 
% z1         = 1000 * 1500;
% z2         = 1588 * 2906;
% [alpha, A] = fx_calcualate_attcoef(Ab, A0, h, z1, z2);
% process3   = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A);
% filename_fig = filename;
% process3.show_surfaces(filename_fig(1:end-5));

% get the surfaces by statistic calculation again
% window       = [1 350 1 80 1 1300];
% [Ab, A0, h]  = process3.statistic_front_backechoes(window);
% z1           = 1000 * 1500;
% z2           = 1588 * 2906;
% [alpha, A]   = fx_calcualate_attcoef(Ab, A0, h, z1, z2);
alpha        = 3e-3;
A_ratio      = 0.5;
% A_ratio      = 0.4; % for 0.3 m drop
process3     = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
process3.show_surfaces(filename_fig(1:end-5));

process3     = process3.recover_surface;
process3     = process3.smooth_rear_I;
process3.show_surfaces(filename_fig(1:end-5));

% tackle the problem of walls determination
TOF_walls    = -mean(process3.front_I-process3.rear_I, 'all', 'omitnan');
% process3.rear_I = process3.front_I + mean(TOF_walls);
% process3.show_surfaces;
% % shift the A scan in the time domain
% process3        = process3.shift_A_scan(200);

%% display A B scans and time gate to determine the back wall

process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, 132, 68);
process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, 150, 150);

show_Raw_filtered_3D(obj, xslice, yslice, zslice)
        
process3.demo_Bscan_interply('x', 20, 1:350, 'healthy');
process3.demo_Bscan_interply('x', 160, 1:350, 'delaminated');

%% original results display
xslice = 178 / process3.fx * 1e3;
yslice = 172 / process3.fy * 1e3;
zslice = []; % us

process3.show_Raw_filtered_3D(xslice, yslice, zslice)
       
process3.show_inaminph_3D(xslice, yslice, zslice);

%% internal damage features
process3.damage_imaging;
%
mid_x     = 178;
mid_y     = 178;
winy      = round(mid_x-23.58e-3*process3.fx/2): round(mid_x+23.58e-3*process3.fx/2);
winx      = round(mid_y-18.67e-3*process3.fy/2): round(mid_y+18.67e-3*process3.fy/2);
% rot_angle = -110;
rot_angle = 20;
process3.show_front_position(winx, winy, rot_angle);

% % 
% centers = [177, 179];
% radii   = 18e-3*process3.fy/2 - 1;
% process3.show_front_position_circle(centers, radii);

%% 3d ply track
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
f0        = 11e6;
sigma     = 0.8;

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
nol      = 49;
% process3 = process3.track_interply_2ndharmonic('img_hil_filter', nol);
process3 = process3.track_interply('img_hil_filter');

% use 2nd-harmonic for first and the last interplies, 
% use fundamental resonance for another interplies
f0_1     = 12.8e6;
sigma0_1 = 0.8;
f0_2     = 6.3e6;
sigma0_2 = 0.7;
nol      = 26; % make it 25+1 in case that there is redundant accidently-tracked interply
process3 = process3.track_interply_hybrid( ...
    'img', f0_1, sigma0_1, f0_2, sigma0_2, nol);

%% display A-scan and hybrid instantaneous phase approach
process3 = process3.track_interply_hybrid_ascan( ...
    'img', f0_1, sigma0_1, f0_2, sigma0_2, nol, 10, 10);

%%
B_type     = 'x';
index      = 300;
Bwin       = 1:100;
TOF_oneply = TOF_walls / 24; % 24 plies
process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

% 3d slices
xslice = 170 / process3.fx * 1e3;
yslice = 190 / process3.fy * 1e3;
xslice = 175 / process3.fx * 1e3;
yslice = 185 / process3.fy * 1e3;
zslice = [];
process3.show_interply_track_3D_knn(xslice, yslice, zslice);

% process.show_track_interply(xslice, yslice, zslice);
win_x      = 1:80;
win_y      = 1:80;
process3.show_oneinterply(05, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(15, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(22, 'logGabor', 3000, win_x, win_y, TOF_oneply);
%
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

%% video of interply tracks
B_type        = 'x';
Bwin          = 1:350;
property_name = 'img_hil';
process3.makeMovie_InterplyTrack(B_type, Bwin, property_name);

%% structrul tensor 
% out-of-plane angles
% logGabor filtered
f0_1     = 6.3e6;
sigma0_1 = 0.7;
process = process3.Filter_logGabor(f0_1, sigma0_1, 'img');

% smoothing scale
ds_rate      = 1;
sigma1       = [3, 3, 3];
% integration scale
sigma2       = [3, 3, 3];
PropertyName = 'img_hil_filter';
process3     = process3.structural_tensor(sigma1, sigma2, PropertyName, ds_rate);

% display the angles 3d slices
medf_kernel  = [3, 3, 3];
xslice       = 220 / process3.fx * 1e3;
yslice       = 170 / process3.fy * 1e3;
zslice       = [];
process3.show_angles_ST_zangle(medf_kernel, xslice, yslice, zslice, ds_rate);

%
B_type     = 'x';
index      = 170;
Bwin       = 1:350;
process3.show_angles_ST_Bscan(B_type, index, Bwin);

%
B_type     = 'y';
index      = 170;
Bwin       = 1:350;
process3.show_angles_ST_Bscan(B_type, index, Bwin);

%% video of out-of-plane angles
B_type        = 'x';
Bwin          = 1:350;
property_name = 'angle_z';
process3.makeMovie_angles_ST_Bscan(B_type, Bwin, property_name);

%% one-plane EDA
% 2D demo
% define the C_scan_inam as well
% % filtered
f0       = 15e6;
sigma    = 0.7; 
process3 = process3.Filter_logGabor(f0, sigma, 'img');

%
PropertyName  = 'img_hil_filter';
% PropertyName = 'img_hil';
% PropertyName = 'img_WienerDeconv';
ply           = 7;
ratio         = 0.5;
process3      = process3.define_plywise_inamCscan(ply, ratio, PropertyName);
% curvelet denoising
sigma         = 5e-3;
%
imagename     = 'C_scan_inam_plywise';
process3      = process3.curvelet_denoise_inamCscan(sigma, imagename);

% RT
theta         = 0:5:180;
radiis        = 16;
% imagename     = 'C_scan_inam_plywise';
% imagename     = 'C_scan_inam_denoise';
process3      = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
angle_compens = 0;
process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);

% 2d log-Gabor fitler
wavelength    = 4:1:16;  
orientation   = 1:1:180; % degree
SFB           = 1;       % [0.5 2.5]
SAR           = 0.5;     % [0.23 0.92]
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K             = 2e0;
process3      = process3.compute_logGabor_filter_withoutFig( ...
    PropertyName, wavelength, orientation, SFB, SAR, imagename);
angle_compens = 20;
process3.show_orientation_by_ID_allwl( ...
    wavelength, orientation, K, imagename, angle_compens);

% % curvelet
% process3  = process3.compute_curvelet(imagename);
% process3.show_orientation_by_curveletID(imagename);

%% multi one-plane EDA 
% 2D demo
theta               = 0:1:180;
ply_wise_Cscan_inam = cell(1, 24);
anguler_1D_ply_wise = cell(1, 24);
center              = [50, 50];
radii               = 49;
PropertyName        = 'img_hil_filter';
imagename           = 'C_scan_inam_plywise';
for ply = 1:24
    ratios  = 0.5:0.1:0.6;
    process3  = process3.define_plywise_inamCscan(ply, ratios, PropertyName);
%     % curvelet denoising
%     sigma_den = 1e-6;
%     process3  = process3.curvelet_denoise_inamCscan(sigma_den, imagename);
    [R, anguler_1D, max_angle_I, xp] = fx_Radonto1Dangular...
        (process3.C_scan_inam_plywise, center, radii, theta);
    % plot the 1D anguler (sum of gradient by central-different kernel)
%     figure, plot(theta, anguler_1D, 'linewidth', 2);
    anguler_1D_ply_wise{ply} = anguler_1D;
    ply_wise_Cscan_inam{ply} = process3.(imagename);
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

%% ply-wise in-plane orientation extraction
% RT
% PropertyName = 'img_hil_filter'; % deconvolved signal
% PropertyName = 'img_WienerDeconv'; % deconvolved signal
PropertyName = 'img_hil_filter'; % orginal signal
radius      = 16;
theta       = 1:1:180; % degree
ratios      = 0.1:0.2:0.9;
plies       = 25;
% sigma_denoi = 5e-4;
sigma_denoi = 0;
process3    = process3.extract_local_orientation_RT_3D_plywise(...
    PropertyName, radius, theta, ratios, plies, sigma_denoi);

% show slice
xslice        = 180 / process3.fx * 1e3;
yslice        = 170 / process3.fy * 1e3;
zslice        = []; % 1.8; % us
angle_compens = 20; % degree
% angle_compens = 0; % degree
process3.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);

process3.show_inplane_direction_3D_meanstd( ...
    xslice, yslice, zslice, plies, ratios, PropertyName, wavelength);

% 2d log-Gabor fitler
wavelength  = 4:2:16;  
orientation = 1:1:180; % degree
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K           = 1e0;
ratios      = 0.1:0.2:0.9;
% K                  = 0;
% for wavelength = [2 4 8 16 32 40]
process3 = process3.extract_local_orientation_3D_plywise_allwl(...
    PropertyName, wavelength, orientation, ratios, plies, K, sigma_denoi);

% show 3d slice
medfiltersize = [3 3 3];
xslice        = 175 / process3.fx * 1e3;
yslice        = 175 / process3.fy * 1e3;
angle_compens = 20; % degree
xslice        = 100 / process3.fx * 1e3;
yslice        = 100 / process3.fy * 1e3;
zslice        = []; % 1.8 us
angle_compens = 0;
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, medfiltersize, angle_compens);

% show bscans in-plane direciton
%
B_type     = 'x';
index      = 178;
Bwin       = 1:350;
process3.show_inplane_direction_3D_ID_Bscan(B_type, index, Bwin, medfiltersize, angle_compens);
%
B_type     = 'y';
index      = 172;
Bwin       = 1:350;
process3.show_inplane_direction_3D_ID_Bscan(B_type, index, Bwin, medfiltersize, angle_compens);

% %
% process3.show_inplane_direction_3D_ID_meanstd( ...
%     xslice, yslice, zslice, plies, ratios, PropertyName, wavelength);
% process3.show_inplane_direction_3D_ID_global(plies, ratios, orientation);

%% video of in-plane angles
B_type        = 'x';
Bwin          = 1:350;
property_name = 'Inplane_direction_3D_ID';
angle_compens = 20; % degree
process3.makeMovie_fiberDirection(B_type, Bwin, property_name, angle_compens);

%% ***********
FolderName = "F:\Xiayang\results\Impact_analysis\";   % the destination folder
% FolderName = "F:\Xiayang\results\Teijin_samples\";
% FolderName = "F:\Xiayang\results\image_processing\";
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '31122020_13h15m49s'));
% '31122020_14h27m51s' Teijin Sample
% '31122020_13h15m49s.mat' % 0.1 m impact sample
% '13072020_16h54m52s' % non-damage sample 
% '14122020_21h50m10s' % 0.3 m damage sample 

%% **************************************
% plus deconvolution
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
% process2.show_hilbert_Ascan(101, 96);
% process2.show_hilbert_Ascan(50, 76);

%% window cutting
% process2.show_img_shape;
window   = [11 110 11 210 100 1400];
window   = [51 350 1 300 1 1400];
x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;
% 
% ratio    = 0.5;
% process3 = process3.time_varying_gain(ratio);

%% read reference signal
[Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
filename               = strcat(Pathname1, Filename1);
x                      = 10;
y                      = 10;
process3               = process3.read_refer(filename, x, y);

% cut the ref. signal
win                    = [1.1e-6 3e-6]; % unit: s
process3               = process3.cut_reference_signal(win);

% calculate the SNR
win_noise              = [0.1e-6 0.5e-6];
win_signal             = [0.6e-6 1e-6]; % unit: s
snr                    = process3.calculate_SNR(win_signal, win_noise);

% align the signal
process3               = process3.align_refer_ascan(x, y); 

%% deconvolution
x           = 100;
y           = 100;
q_factor    = 1e-2;
q_factor_AR = 1e-2;
fft_padding = 2^11;
f1          = 8e6;
f2          = 19e6;
bw          = [f1 f2];
bwr         = -3:-1:-10;
k           = 45;
lam         = 0.01;
ker_wid     = 200;
Nit         = 45;
DownRate    = 1;
fig_subname = '2-CEH105-24-p4_15MHz';
process3.demo_deconvolutions(x, y, q_factor, q_factor_AR, ...
    fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit, fig_subname);

% demo B scan
B_type      = 'x'; % scaning direction 'y'
index       = 50;
Bwin        = 1:100;
% process = process.deconv_ARextrapolation(Bwin, B_type, index, q_factor_AR, fft_padding, bw, bwr, k, ker_wid, fig_subname);

% process.demo_deconvolutions(x, y, q_factor, q_factor_AR, fft_padding, bw, k, lam, ker_wid, DownRate, Nit, fig_subname);
% 
process3     = process3.apply_deconvolutions_Bscan(Bwin, B_type, index, q_factor, q_factor_AR, fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit);
process3.demo_deconvolutions_Bscans(B_type, index, Bwin, fig_subname);

%%
process3      = process3.apply_deconvolutions_onlyWiener(q_factor, ker_wid);

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
MinPD        = 500;
MinPH        = 0.03;  % these 2 parameters need to be changed for surface estimation.
% PropertyName = 'img_WienerDeconv'; bh
PropertyName = 'img_hil';
process3.show_Ascan_inam_peaks(128, 148, MinPD, MinPH, PropertyName); % x, y
% surface calculation
max_len                = 900;
threshold_delamination = 0.5;
process3               = process3.find_front_amp(...
    MinPD, MinPH, PropertyName, max_len, threshold_delamination);
filename_fig           = filename;
process3.show_surfaces(filename_fig(1:end-5));
%
process3.smooth_rear_I;
process3.show_surfaces(filename_fig(1:end-5));
% %
% load TOF.mat;
% process2.rear_I = process2.front_I + mean(TOF_walls);
% % tackle the problem of walls determination
% process2.show_surfaces;

%% original results display
xslice = 10 / process3.fx * 1e3;
yslice = 10 / process3.fy * 1e3;
zslice = 2.4; % us
process3.show_inaminph_3D(xslice, yslice, zslice);

%% ply track by inam (deconvolve)
property_name = 'img_WienerDeconv';
% property_name = 'img_hil';
process3      = process3.track_interply_inam(property_name, 25);

B_type     = 'x';
index      = 10;
Bwin       = 1:100;
process3.show_B_scan_interply(B_type, index, Bwin, property_name);

% show ply track
win_x  = 1:100;
win_y  = 1:100;
TOF_oneply = TOF_walls / 24; % 24 plies
process3.show_oneinterply(05, property_name, 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(15, property_name, 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(22, property_name, 3000, win_x, win_y, TOF_oneply);

process3.thickness_estimation(5.5, 24, win_x, win_y, std_mean_threshold, '15MHz');

% plot with uniformly spaced locaitons
nol                = 25;
TOF                = process3.track_show_B_scan_interply(...
    B_type, index, Bwin, property_name, threshold, nol, std_mean_threshold);

%% 3d ply track
% add noise to img.
% process = process.addnoise(25);
% origin without filter
threshold  = 0.01;  % max(inam) * threshold for distinguishing the back surface

% % filtered
f0        = 10.89e6;
sigma     = 0.8;
% for debug
% threshold = 0.05;
x         = 30;
y         = 30;
% process3.demo_logGabor_plytrack_inph_v2(x, y, f0, sigma, threshold);
process3.show_hilbert_Ascan(x, y);
process3.show_logGabor_Scaleogram(5e6:0.1e6:15e6, sigma, x, y);
%  
process3 = process3.Filter_logGabor(f0, sigma, 'img');
process3 = process3.track_interply('img_hil_filter');
% 3d slices
xslice   = 125 / process3.fx * 1e3;
yslice   = 125 / process3.fy * 1e3;
zslice   = [];
process3.show_interply_track_3D_knn(xslice, yslice, zslice);
% 
B_type     = 'y';
index      = 125;
Bwin       = 1:250;
TOF_oneply = TOF_walls / 24; % 24 plies
process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil_filter');

%% structrul tensor 
% out-of-plane angles
% logGabor filtered
f0       = 10.9e6;
sigma    = 0.8;
process3 = process3.Filter_logGabor(f0, sigma, 'img');

% smoothing scale
ds_rate      = 1;
sigma1       = [3, 3, 3];
% integration scale
sigma2       = [3, 3, 3];
PropertyName = 'img_hil_filter';
process3     = process3.structural_tensor(sigma1, sigma2, PropertyName, ds_rate);

% % display the angles 3d slices
% medf_kernel  = [3, 3, 3];
% xslice       = 175 / process3.fx * 1e3;
% yslice       = 175 / process3.fx * 1e3;
% zslice       = [];
% process3.show_angles_ST_zangle(medf_kernel, xslice, yslice, zslice, ds_rate);

%
B_type     = 'x';
index      = 125;
Bwin       = 1:250;
process.show_angles_ST_Bscan(B_type, index, Bwin);

%
B_type     = 'y';
index      = 125;
Bwin       = 1:250;
process.show_angles_ST_Bscan(B_type, index, Bwin);

%% one-plane EDA
% 2D demo
% total_plies = 24;
% [angle, process] = process.show_one_Radontransform_vitualprofile(5, 0.5, 'img_hil_filter', center, radii, total_plies);

% % filtered
f0       = 15e6;
sigma    = 0.7; 
process3 = process3.Filter_logGabor(f0, sigma, 'img');

% define the C_scan_inam as well
% PropertyName = 'img_WienerDeconv';
% PropertyName = 'img_hil_filter';
PropertyName = 'img_hil';
ratio        = 23.5/24;
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
angle_compens = -7;
process3      = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);
% 2d log-Gabor fitler
wavelength  = 4:2:20;
orientation = 1:1:180;
SFB         = 1; % [0.5 2.5]
SAR         = 0.7; % [0.23 0.92]
% imagename = 'C_scan_inam_denoise';
% process3     = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% process3.show_orientation_by_ID_allwl(wavelength, orientation, K);
imagename     = 'C_scan_inam';
angle_compens = -7;
process3      = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K = 0.5e0;
process3.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);

%% 'distance to front and rear' in-plane orientation extraction
% RT
% PropertyName     = 'img_WienerDeconv'; % deconvolved signal
PropertyName  = 'img_hil_filter'; % orginal signal
% PropertyName  = 'img_hil';

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
ratio              = 0/24:0.2/24:24/24;
K                  = 1e0;
% sigma              = 5e-4;
sigma              = 0;
tic;
process3           = process3.extract_local_orientation_3D_parallel_allwl(...
    PropertyName, wavelength, orientation, ratio, K, sigma);
toc;

% show slice
xslice        = 150 / process3.fx * 1e3;
yslice        = 50 / process3.fy * 1e3;
zslice        = [];
mfsize        = [1 1 1];
angle_compens = -7;
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

% 
B_type = 'x';
index  = 150;
Bwin   = 1:300;
process3.show_inplane_direction_3D_ID_Bscan(B_type, index, Bwin, [3 3 3], angle_compens);
%
B_type = 'y';
index  = 50;
Bwin   = 1:300;
process3.show_inplane_direction_3D_ID_Bscan(B_type, index, Bwin, [3 3 3], angle_compens);

%% ************** calculate the mean fiber angle and its standard deviation *************
% need the reference angle to calcualate the mean and std!
Stacking_sequence = [
    45 0 -45 -90 ...
    45 0 -45 -90 ...
    45 0 -45 -90 ...
    -90 -45 0 45 ...
    -90 -45 0 45 ...
    -90 -45 0 45];
m_fiber_angle_arr   = NaN(1, 24);
std_fiber_angle_arr = NaN(1, 24);
Idof_N_arr          = NaN(24, 180);
PropertyName_IFD    = 'Inplane_direction_3D_ID';
for p = 1:24 % 24 layers for default
    if p ==24
        ratio = (p-0.5) * 1/24;
    else
        ratio = (p-0.5) * 1/24;
    end
    ratio_next = (p-0.3) * 1/24;
    %     ref_angle  = Stacking_sequence(p);
    Idof_N_arr(p,:) = process3.calculate_angle_distrib(PropertyName, PropertyName_IFD, ratio, angle_compens);
    %     std_fiber_angle_arr(p) = std_fiber_angle_arr(p) / sum(Idof_Nr_arr(p,:));
    % normalize
    %     Idof_N_arr(p,:) = (Idof_N_arr(p,:)) / sum(Idof_N_arr(p,:));
    %     % log
    %     Idof_N_arr(p,:) = log(Idof_N_arr(p,:));
    disp([num2str(p) '/' '24']);
    if p==3 || p==11 || p==22
        figure('Name', ['angle_distribution_' num2str(p) '_' PropertyName_IFD]);
        set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
        h1 = bar(-89-22:90-22, Idof_N_arr(p,:)/max(Idof_N_arr(p,:)));
        hold on;
        xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
        ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
        %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
        set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
        set(gca, 'linewidth', 2);
%         set(gca, 'YScale', 'log')
        %         legend([h1 h2], 'Angle distribution', 'Multiple Gaussian fitting', 'Interply track');
        %         ytickformat('percentage');
    end
end
%
maps             = cell(1,24);
for i = 1:24
    maps{i} = ['Ply' num2str(i) ':    ' num2str(Stacking_sequence(i)) '\circ'];
end
fx_creat_chartPlot(maps, Idof_N_arr, Stacking_sequence);

% m_fiber_angle_arr = m_fiber_angle_arr-90-23+angle_compens;
% 
m_fiber_angle_arr_RT   = NaN(1, 24);
std_fiber_angle_arr_RT = NaN(1, 24);
Idof_N_arr             = NaN(24, 180);
PropertyName_IFD       = 'Inplane_direction_3D';
angle_compens          = -7;
for p = 1:24 % 24 layers for default
    if p ==24
        ratio = (p-0.5) * 1/24;
    else
        ratio = (p-0.5) * 1/24;
    end
    ratio_next = (p-0.3) * 1/24;
    ref_angle  = Stacking_sequence(p);
    %     [m_fiber_angle_arr_RT(p), std_fiber_angle_arr_RT(p), yhat, Idof_N_arr(p,:)] = process3.calculate_m_std_fiber_angle( ...
    %         PropertyName, PropertyName_IFD, ratio, ratio_next, angle_compens);
    %     std_fiber_angle_arr_RT(p) = std_fiber_angle_arr_RT(p) / sum(Idof_N_arr(p,:));
    Idof_N_arr(p,:) = process3.calculate_angle_distrib(PropertyName, PropertyName_IFD, ratio, angle_compens);
    %
    %     % normalize
    %     Idof_N_arr(p,:) = (Idof_N_arr(p,:)) / sum(Idof_N_arr(p,:));
    %     % log
    %     Idof_N_arr(p,:) = log(Idof_N_arr(p,:));
    disp([num2str(p) '/' '24']);
    if p==3 || p==11 || p==22
        figure('Name', ['angle_distribution_' num2str(p) '_' PropertyName_IFD]);
        set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
        h1 = bar(-89-22:90-22, Idof_N_arr(p,:)/max(Idof_N_arr(p,:)));
        xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
        ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
        %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
        set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
        set(gca, 'linewidth', 2);
        %         set(gca, 'YScale', 'log')
        %         legend([h1 h2], 'Angle distribution', 'Multiple Gaussian fitting', 'Interply track');
        %         ytickformat('percentage');
    end
end
fx_creat_chartPlot(maps, Idof_N_arr, Stacking_sequence);

% m_fiber_angle_arr_RT = m_fiber_angle_arr_RT-90-23 + angle_compens;
% 
% T = array2table([m_fiber_angle_arr; m_fiber_angle_arr_RT; ...
%     std_fiber_angle_arr; std_fiber_angle_arr_RT]);
% writetable(T,'test.xlsx','Sheet',1);

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

% saveas(gcf,'Barchart','epsc')

%% save all figures 4th paper
FolderName  = "F:\Xiayang\results\Impact_analysis";   % the destination folder
% FolderName  = "F:\Xiayang\results\Impact_analysis\ABCscans";   % the destination folder
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

%% save to bin file

% matrix = process3.img;
% matrix = single(matrix);
% save('output.mat', 'matrix');
% 
% load('output.mat');
