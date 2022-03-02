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

process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, 10, 10);
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

%% multi slices EDA

PropertyName  = 'Inplane_direction_3D_ID';
% PropertyName  = 'angle_z';
cf1 = figure('Name', ['Cscan_slice_' , PropertyName]);
set(cf1, 'Position', [0, 0, 600, 600], 'color', 'white');
for ratio = [6.5 13.5 20.5]/24
% for ratio = [3.5 4.5 9.5 14.5]/24
    [process3, ~, C_scan_index]    = process3.define_parallel_inamCscan...
        (ratio, PropertyName);
    inam = process3.display_slice(C_scan_index, PropertyName);
    [lx, ly] = size(C_scan_index);
    X = (1: lx) / process3.fx * 1e3;
    Y = (1: ly) / process3.fy * 1e3;
    C_scan_index = C_scan_index / process3.fs * 1e3 * 3000 / 2;
    figure(cf1)
    surf(X, Y, C_scan_index, inam, 'EdgeColor', 'none');
    set(gca, 'Zdir', 'reverse');
    hold on;
    % angle
    colormap(hsv);
%     caxis([0 2]);
%     % direction
%     colormap(hsv);
%     caxis([-90 90]);
    colorbar;
    xlim([0 X(end)]);
    ylim([0 Y(end)]);
    xlabel('\fontname {times new roman} X (mm)', 'fontsize', 16);
    ylabel('\fontname {times new roman} Y (mm)', 'fontsize', 16);
    zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
end

%%  thickness multi slices
cf1 = figure('Name', ['Cscan_slice_' , 'thickness']);
set(cf1, 'Position', [0, 0, 600, 600], 'color', 'white');
for ratio = [6.5 13.5 20.5]/24
    [~, ~, C_scan_index]    = process3.define_parallel_inamCscan...
        (ratio, PropertyName);
    C_scan_index = C_scan_index / process3.fs * 1e3 * 3000 / 2;
    inam = process3.est_thickness(:, :, ceil(ratio*24));
    inam = inam / process3.fs * 1e6 * 3000 / 2;
    [lx, ly] = size(inam);
    X = (1: lx) / process3.fx * 1e3;
    Y = (1: ly) / process3.fy * 1e3;
    figure(cf1)
    surf(X, Y, C_scan_index, inam, 'EdgeColor', 'none');
    set(gca, 'Zdir', 'reverse');
    hold on;
    colormap(jet);
    colorbar;
    caxis([100 300]);
    xlim([0 X(end)]);
    ylim([0 Y(end)]);
    xlabel('\fontname {times new roman} X (mm)', 'fontsize', 16);
    ylabel('\fontname {times new roman} Y (mm)', 'fontsize', 16);
    zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
end

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

%% structrul tensor 
% out-of-plane angles
% logGabor filtered
f0_1     = 12.8e6;
sigma0_1 = 0.8;
process3 = process3.Filter_logGabor(f0_1, sigma0_1, 'img');

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

%% one-plane EDA

% define the C_scan_inam as well
% PropertyName = 'img_WienerDeconv';
% PropertyName = 'img_hil_filter';
PropertyName  = 'img_hil';
angle_compens = -7;
for ratio = [6.5 13.5 20.5]/24
    % 2d log-Gabor fitler
    wavelength  = 4:2:20;
    orientation = 1:1:180;
    SFB         = 1; % [0.5 2.5]
    SAR         = 0.7; % [0.23 0.92]
    % imagename = 'C_scan_inam_denoise';
    % process3     = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
    % process3.show_orientation_by_ID_allwl(wavelength, orientation, K);
    process3      = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
    % K: controls how much smoothing is applied to the Gabor magnitude responses.
    K = 0.5e0;
    process3.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);
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

%% ***********
FolderName = "F:\Xiayang\results\Impact_analysis\";   % the destination folder
% FolderName = "F:\Xiayang\results\Teijin_samples\";
% FolderName = "F:\Xiayang\results\image_processing\";
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '13072020_16h54m52s.mat'));
% '31122020_14h27m51s' Teijin Sample
% '31122020_13h15m49s.mat' % impact sample
% '13072020_16h54m52s' % non-damage sample 


%% save all figures
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

