% read the file

clc;
close all;
fclose all;
clear;

%% Define the class
% read the preprocessed .mat
[Filename1, Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename               = strcat(Pathname1, Filename1);
process3               = class_process_woven_RFdata(filename);

% input the settings 
process3.fs = 250e6; % MHz
process3.fx = 1 / 0.2e-3; % unit: 1 / m
process3.fy = 1 / 0.2e-3; % unit: 1 / m
% process3.fx = 1 / 0.25e-3; % unit: 1 / m
% process3.fy = 1 / 0.25e-3; % unit: 1 / m

% % read the settings from excel
% [Filename1, Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
% settinig_file          = strcat(Pathname1, Filename1);
% 
% % read the data from the struct, or reset the dataset
% process3               = process3.loadSettings(settinig_file);

process3 = process3.read_origin_data; % read (reset) the dataset

%% ***********
FolderName = "F:\Xiayang\results\Woven_samples\";   % the destination folder
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '20220428_BP1_3_025m_V313_25db_PEmat.mat'));
% '20220207_BP1_3_025m_V313_25db_PEmat.mat'
% '20220303_BP1_3_025m_h5m_10db_PEmat.mat'
% '20220428_BP1_3_025m_V313_25db_filp_PEmat.mat'
% '20220428_BP1_3_025m_V313_25db_PEmat'
% '20220428_BP1_3_025m_V313_47db_PEmat'

%% EDA
% orthosliceViewer(process3.img);

process3.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500]; 
% window   = [1 296 1 175 1 1200]; 
window   = [171 270 1 100 1 1500];  % for 025m impact samples
% window   = [61 140 57 136 1 1500];  % for 025m impact samples

x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;
% % 3D viewer
% volumeViewer(abs(process3.img_hil));

%% show A scan
% define the index to select Ascan 
x   = 56;
y   = 63;
% process3.demo_AS_3D_inclinfq(x, y);
process3.show_hilbert_Ascan(x, y);

%% show C scan
% define the index to select Ascan 
close all;
for z = 1150:50:1150
    % z = ;
    PropertyName = 'img';
    process3 = process3.show_Cscan(z, PropertyName);
%     process3.check2dfft_inclphase(z, PropertyName); % 2d spectrum
end

%% 2nd fft viewer
temp = process3.img;
temp = fft(temp, [], 3);
temp = abs(temp(:, :, 1:end/2));
temp = fft(temp, [], 3);
temp = abs(temp(:,:,1:end/2));
fx_Scrollable_3d_view(temp);

%% low-pass filter
bandpassFreq = 30e6;
process3     = process3.Filter_lowpass(bandpassFreq, 'img_hil');

%% 3D volume filter
sigma    = [0.5 0.5 1.5];
process3 = process3.Filter_3DGaussian(sigma, 'img_hil');

%% surface search
global min_pks 
min_pks      = 0.4;
PropertyName = 'img_hil';
% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y
clc;
% close all;

MinPD   = 5;
MinPH   = 0.12;  % these 2 parameters need to be changed for surface estimation.
% surface calculation
max_len = 1500;
alpha   = 4e-3;
A_ratio = 0.95;

process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, ...
    23, 35);

%%
% A_ratio      = 0.4; % for 0.3 m drop
process3     = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
close all;

process3.smooth_rear_I;
process3.show_surfaces(filename_fig(1:end-5));

process3 = process3.recover_surface;

%% normal time window
PropertyName = 'img_hil';
% PropertyName = 'img_hil_filter';
delay        = 600;
max_len      = 1500;
flag         = 0;
front_I_max  = 400;
process3     = process3.find_damage_timewin_asignsurface(PropertyName, max_len, flag, delay, front_I_max, 1);
process3.show_surfaces(filename_fig(1:end-5));

%% align the front surfaces
delay_I      = 300;
PropertyName = 'img_hil';
% PropertyName = 'img';
process3     = process3.align_front_surface(delay_I, PropertyName);

% back surface
% process3     = process3.align_rear_surface(delay_I, PropertyName);

%%
temp = process3.fft2d_mask_fingerprint;

temp_abs = squeeze(temp(:,:,:,1));
% fill nan  
temp_abs(isnan(temp_abs)) = min(temp_abs(:));
orthosliceViewer(temp_abs, 'ScaleFactors', [1 1 1], 'Colormap', jet);

temp_abs = squeeze(temp(:,:,:,2));
% fill nan  
temp_abs(isnan(temp_abs)) = min(temp_abs(:));
orthosliceViewer(temp_abs, 'ScaleFactors', [1 1 1], 'Colormap', jet);

%process3.fft2d_abs_fingerprint;
close all;
figure, subplot(2, 2, 1);
imagesc(squeeze(sum(temp_abs, 3, 'omitnan')));
subplot(2, 2, 2);
imagesc(squeeze(sum(temp_abs, 1, 'omitnan')));
subplot(2, 2, 3);
temp_yz = squeeze(sum(temp_abs, 2, 'omitnan'));
imagesc(temp_yz.');

temp_pha = temp(:,:,:,2);
% process3.fft2d_pha_fingerprint;
figure, subplot(2, 2, 1);
imagesc(squeeze(sum(temp_pha, 3, 'omitnan')));
subplot(2, 2, 2);
imagesc(squeeze(sum(temp_pha, 1, 'omitnan')));
subplot(2, 2, 3);
temp_yz = squeeze(sum(temp_pha, 2, 'omitnan'));
imagesc(temp_yz.');

% imagesc(temp(:,:, 50));
%% demonstrate A-scan and ply track by inph 
% logGabor fitler added
threshold = 0.05;
f0        = 6.3e6;
sigma     = 0.7;
x = 10;
y = 10;
process3.demo_logGabor_plytrack_inph(x, y, f0, sigma, threshold);

%% 3d ply track

% log-Gabor filter
f0        = 12.8e6;
sigma     = 0.7;
process3  = process3.Filter_logGabor(f0, sigma, 'img_hil');

% % low-pass filter
% bandpassFreq = 10e6;
% process3     = process3.Filter_lowpass(bandpassFreq, 'img_hil');

threshold = 0.05;
% process3  = process3.track_interply('img_hil_filter');
process3  = process3.track_interply_inph(threshold, 'img_hil_filter', 21);

%%
% use 2nd-harmonic for first and the last interplies, 
% use fundamental resonance for another interplies
f0_1     = 12.8e6;
sigma0_1 = 0.8;
f0_2     = 6.3e6;
sigma0_2 = 0.7;
nol      = 21; % make it 25+1 in case that there is redundant accidently-tracked interply
process3 = process3.track_interply_hybrid( ...
    'img', f0_1, sigma0_1, f0_2, sigma0_2, nol);

%%
TOF_walls    = -mean(process3.front_I - process3.rear_I, 'all', 'omitnan');

close all;
B_type     = 'x';
index      = 40;
Bwin       = 1:80;
TOF_oneply = TOF_walls / 20; % 24 plies
process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

% % process.show_track_interply(xslice, yslice, zslice);
% win_x = 1:270;
% win_y = 1:100;
process3.show_oneinterply(3, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply(10, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply(15, 'logGabor', 3000, win_x, win_y, TOF_oneply);

%%
% close all;
win_x = 1:100;
win_y = 1:100;

process3.show_oneinterply_2dfft(7, 'logGabor', 3000, win_x, win_y, TOF_oneply);

process3.show_oneinterply_2dfft(13, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply_2dfft(14, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply_2dfft(15, 'logGabor', 3000, win_x, win_y, TOF_oneply);

close all;
process3.show_oneinterply_2dfft(9, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply_2dfft(10, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply_2dfft(11, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3 = process3.show_oneinterply_2dfft(14, 'logGabor', 3000, win_x, win_y, TOF_oneply);

%% save images 
PropertyName = 'img_hil';
ratios       = 0.2:0.1:0.8;
sigma        = 5e-3;
imagename    = 'C_scan_inam_plywise';

[process3, image1, ~]  = process3.define_plywise_inamCscan(1, ratios, PropertyName);
process3.curvelet_denoise_inamCscan(sigma, imagename);
[process3, image2, ~]  = process3.define_plywise_inamCscan(2, ratios, PropertyName);
process3.curvelet_denoise_inamCscan(sigma, imagename);
[process3, image3, ~]  = process3.define_plywise_inamCscan(3, ratios, PropertyName);
process3.curvelet_denoise_inamCscan(sigma, imagename);
[process3, image4, ~]  = process3.define_plywise_inamCscan(4, ratios, PropertyName);
process3.curvelet_denoise_inamCscan(sigma, imagename);

% save('wovenstrcuture_5MHz.mat','image1','image2','image3','image4');

%% save images 
imagename = 'C_scan_inam';

[process3, image1, ~]  = process3.define_parallel_inamCscan(1/20, PropertyName);
process3.curvelet_denoise_inamCscan(sigma, imagename);
[process3, image2, ~]  = process3.define_parallel_inamCscan(2/20, PropertyName);
process3.curvelet_denoise_inamCscan(sigma, imagename);
[process3, image3, ~]  = process3.define_parallel_inamCscan(3/20, PropertyName);
process3.curvelet_denoise_inamCscan(sigma, imagename);
[process3, image4, ~]  = process3.define_parallel_inamCscan(4/20, PropertyName);
process3.curvelet_denoise_inamCscan(sigma, imagename);

save('wovenstrcuture_5MHz.mat','image1','image2','image3','image4');

%% read reference signal
[Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
filename               = strcat(Pathname1, Filename1);
x                      = 10;
y                      = 10;
process3               = process3.read_refer(filename, x, y);

% cut the ref. signal
win                    = [0.1e-6 2e-6]; % unit: s
process3               = process3.cut_reference_signal(win);

% calculate the SNR
win_noise              = [0.1e-6 0.5e-6];
win_signal             = [0.6e-6 1e-6]; % unit: s
snr                    = process3.calculate_SNR(win_signal, win_noise);

% align the signal
process3               = process3.align_refer_ascan(x, y); 

%% deconvolution
x           = 50;
y           = 50;
q_factor    = 1e-2;
q_factor_AR = 1e-2;
fft_padding = 2^11;
f1          = 8e6;
f2          = 19e6;
bw          = [f1 f2];
bwr         = -3:-1:-10;
k           = 30;
lam         = 0.01;
ker_wid     = 100;
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
% process3 = process3.apply_deconvolutions(q_factor, q_factor_AR, fft_padding, bw, bwr, k, ker_wid);
process3 = process3.apply_deconvolutions_onlyWiener(q_factor, ker_wid);

% save('process2.mat','process2');

fx_Scrollable_3d_view(abs(process3.img_hil));

% fx_Scrollable_3d_view(abs(fft(abs(fft(process3.img, [], 3)), [],3)));

%% one-plane EDA
% % 2D demo
% % define the C_scan_inam as well
% % filtered
% f0       = 15e6;
% sigma    = 0.7; 
% process3 = process3.Filter_logGabor(f0, sigma, 'img');
% low-pass filter
% bandpassFreq = 30e6;
% process3     = process3.Filter_lowpass(bandpassFreq, 'img_hil');

%
% PropertyName = 'img_hil_filter';
PropertyName = 'img_hil';
% PropertyName = 'img_WienerDeconv';
% PropertyName = 'img_WienerDeconv_AR';
imagename   = 'C_scan_inam';

% ply    = 7;
% ratios = [0.2:0.1:0.8];
% ratios = 0.1;
% process3  = process3.define_plywise_inamCscan(ply, ratios, PropertyName);
% imagename = 'C_scan_inam_plywise';
process3 = process3.define_parallel_inamCscan(3.5/20, PropertyName);
% curvelet denoising
sigma    = 5e-3;
%
process3  = process3.curvelet_denoise_inamCscan(sigma, imagename);

% % ********** for EDA
% for ply = 1:1
%     process3 = process3.define_plywise_inamCscan(ply, ratio, PropertyName);
%     imagename = 'C_scan_inam_plywise';
%     process3  = process3.curvelet_denoise_inamCscan(sigma, imagename);
% end
% % *************

% % RT
% theta         = 1:1:180;
% radiis        = 10;
% % imagename     = 'C_scan_inam_plywise';
% % imagename     = 'C_scan_inam_denoise';
% process3      = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
% angle_compens = 0;
% process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);

% process3 = process3.show_oneinterply_2dfft(5, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% imagename   = 'mono_img';
% 2d log-Gabor fitler
wavelength  = 8:2:12;
orientation = 0:2:180; % degree
SFB         = 1.5;       % [0.5 2.5]
SAR         = 0.5;     % [0.23 0.92]
process3    = process3.compute_logGabor_filter_withoutFig( ...
    PropertyName, wavelength, orientation, SFB, SAR, imagename);
angle_compens = 0;
% close all
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K             = 1e-1;
process3.show_orientation_by_ID_allwl( ...
    wavelength, orientation, K, imagename, angle_compens);

% % curvelet
% process3  = process3.compute_curvelet(imagename);
% process3.show_orientation_by_curveletID(imagename);

%% monogenic signal analysis
PropertyName = 'img_hil';
process3     = process3.define_parallel_inamCscan(3.5/20, PropertyName);
imagename    = 'C_scan_inam';
% fillna
C_scan_inam = process3.(imagename);
C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');

% ******** get the image to negative and positive periodic
C_scan_inam = C_scan_inam - mean(C_scan_inam);

%
[Y, X]       = size(C_scan_inam);
cw           = 10;
filtStruct   = createMonogenicFilters(Y, X, cw, 'lg', 0.55);
[m1, m2, m3] = monogenicSignal(C_scan_inam, filtStruct);
% [m1, m2, m3] = monogenicSignal_nofilter(C_scan_inam, filtStruct);

% close all;
% Local energy (calculated on a per-scale basis)
LE = localEnergy(m1, m2, m3);
% Local phase (calculated on a per-scale basis)
LP = localPhase(m1, m2, m3);
% Local orientation (calculated on a per-scale basis)
% Only need to pass the odd parts (m2,m3) as even part (m1) is irrelevant
LO = localOrientation(m2, m3);
%
close all;
% Display
figure;
imagesc(C_scan_inam), axis image; colormap(jet);
h = colorbar;
set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
% colormap gray
title('Test Image');

figure();
imagesc(reshape(LE, Y, [])), axis image;
h = colorbar;
colormap(jet);
set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
title('Local Energy Over Scales');
figure();
imagesc(reshape(LP, Y, [])), axis image, axis off, colormap gray;
title('Local Phase Over Scales');
figure();
imagesc(reshape(LO, Y, [])), axis image, axis off; colorbar;
title('Local Orientation Over Scales (radians)');
h = colorbar;
colormap(hsv);

%% 3D monogenic signal analysis
clc;
D         = process3.img;
[Y, X, Z] = size(D);

% woven strucutre
cw_xy = 10;
cw_z  = 250e6/12.8e6;
% D     = imgaussfilt3(D, 1);
% Now use these wavelengths to create a structure containing
% frequency-domain filters to calculate the monogenic signal. Same as for
% the 2D case, except now we need to pass three image dimensions and the
% two sets of wavelengths
filtStruct = createMonogenicFilters3D(Y, X, Z, cw_xy, cw_z, 'lg', 0.55);
% Now we can use this structure to find the monogenic signal for the volume
[m1,m2,m3,m4] = monogenicSignal3D(D, filtStruct);
% Local energy (calculated on a per-scale basis)
LE = localEnergy3D(m1,m2,m3,m4);
% Local amplitude (calculated on a per-scale basis)
LA = sqrt(LE);
% Local phase (calculated on a per-scale basis)
LP = localPhase3D(m1,m2,m3,m4);
[FS,FA] = featureSymmetry3D(m1,m2,m3,m4);

% PC = phaseCongruency3D(m1,m2,m3,m4,0.05);

% % Cscan image
% z  = 400;
% im = LA(:, :, z);
% imagesc(im), axis image; colormap(jet);
% h = colorbar;
% set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
% 
% % Ascan
% x = 45;
% y = 45;
% figure,plot(squeeze(LA(x, y, :)));

LP_1D = abs(process3.img_hil);

fx_Scrollable_3d_view(LP_1D);

%% multi one-plane EDA 
% 2D demo
ply_num             = 20;
ply_wise_Cscan_inam = cell(1, ply_num);
Idof_N_arr          = NaN(ply_num, 180);
% PropertyName        = 'img_hil';
% PropertyName        = 'img_hil_filter';
PropertyName        = 'img_WienerDeconv';
% imagename           = 'C_scan_inam_plywise';
imagename           = 'C_scan_inam';

for ply = 1:ply_num
%     process3     = process3.show_oneinterply_2dfft(ply, 'logGabor', 3000, win_x, win_y, TOF_oneply);
%     process3 = process3.define_plywise_inamCscan(ply, ratios, PropertyName);
    process3 = process3.define_parallel_inamCscan(ply/2/ply_num, PropertyName);
    process3 = process3.compute_logGabor_filter_withoutFig( ...
        PropertyName, wavelength, orientation, SFB, SAR, imagename);
    [image_ori, Idof_N] = process3.show_orientation_by_ID_allwl( ...
        wavelength, orientation, K, imagename, angle_compens);
    ply_wise_Cscan_inam{ply} = image_ori;
    Idof_N_arr(ply, :) = Idof_N;
    close all;
end

% display all images in on figure
close all;
figure('Name', 'ply_wise_Cscan_inam')
for ply = 1:ply_num
    subplot(ply_num/10, 10, ply)
    pcolor(ply_wise_Cscan_inam{ply}); shading flat;
    colormap hsv;
    caxis([-90 90]);
    title(num2str(ply));
end
colorbar;

% **** 
Stacking_sequence = [
    45 0 -45 -90 ...
    45 0 -45 -90 ...
    45 0 0 45 ...
    -90 -45 0 45 ...
    -90 -45 0 45];
maps = cell(1, ply_num);
for i = 1: ply_num
    maps{i} = ['Ply' num2str(i) ':    ' num2str(Stacking_sequence(i)) '\circ'];
end
%
Idof_N_arr = circshift(Idof_N_arr, 22, 2);
fx_creat_chartPlot(maps, Idof_N_arr, Stacking_sequence);
%*********

%
figure('Name', 'ply_wise_anguler_1D')
anguler_stacking = nan(180, ply_num);
for ply = 1:ply_num
    subplot(ply_num/4, 4, ply)
    plot(-89:90, anguler_1D_ply_wise{ply}, 'linewidth', 2);
    anguler_stacking(:, ply) = anguler_1D_ply_wise{ply};
end

fx_Scrollable_3d_view(abs(process3.img_hil));

%% 'distance to front and rear' in-plane orientation extraction
% RT
PropertyName  = 'img_hil';
% PropertyName  = 'img_hil_filter';
% PropertyName = 'img_WienerDeconv';

% radius        = 10;
% theta         = 1:1:180; % degree
ratio         = 0/20:0.1/20:20/20;
% sigma_denoise = 0;
% tic;
% process3      = process3.extract_local_orientation_RT_3D_parallel(...
%     PropertyName, radius, theta, ratio, sigma_denoise);
% toc;
% 
% % show slice
% xslice        = 100 / process3.fx * 1e3;
% yslice        = 100 / process3.fy * 1e3;
% zslice        = []; % us
% angle_compens = 0;
% process3.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);

% 2D log-Gabor filter
wavelength  = 10:2:10;
orientation = 0:1:180; % degree
SFB         = 1.5;       % [0.5 2.5]
SAR         = 0.5;     % [0.23 0.92]
K           = 0;   
% sigma              = 5e-4;
sigma       = 0;
tic;
process3 = process3.extract_local_orientation_3D_parallel_allwl(...
    PropertyName, wavelength, orientation, ratio, K, sigma);
toc;


% show slice
xslice        = 50 / process3.fx * 1e3;
yslice        = 50 / process3.fy * 1e3;
zslice        = [];
mfsize        = [1 1 1];
angle_compens = -14;
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

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

%% trial on cepstrum
temp_img  = process3.img;
[m, n, l] = size(temp_img);
temp_img_ceps = nan(m, n, l);
for i = 1:m
    for j = 1:n
        ascan            = squeeze(temp_img(i, j, :));
        ascan_fft        = fft(ascan, 2*l);
        ascan_fft        = ascan_fft(1:end/2);
        ascan_fft        = log(abs(ascan_fft)) + 1i*angle(ascan_fft);
        ascan_fft2       = ifft(ascan_fft, 2*l);
        ascan_fft2(1:15) = 0;
        % inverse
        ascan_fft = fft(ascan_fft2(1:end/2), 2*l);
        ascan_fft = ascan_fft(1: end/2);
        ascan_fft = exp(ascan_fft);
        ascan     = ifft(ascan_fft, l);
        temp_img_ceps(i, j, 1:length(ascan)) = ascan;
    end
end
    
fx_Scrollable_3d_view(real(temp_img_ceps));


%%
D_3d = fx_Scrollable_3d_view(process3.Inplane_direction_3D_ID);

process3.statistic_angular_distribution(angle_compens);

%% video of in-plane angles
B_type        = 'x';
Bwin          = 1:100;
property_name = 'Inplane_direction_3D_ID';
process3.makeMovie_fiberDirection(B_type, Bwin, property_name, angle_compens);

%% ply-wise in-plane orientation extraction
% RT
% PropertyName = 'img_hil_filter'; % deconvolved signal
% PropertyName = 'img_WienerDeconv'; % deconvolved signal
PropertyName = 'img_hil'; % orginal signal
wavelength   = 8:2:28;
orientation  = 0:10:180;
ratio        = 0/20:0.1/20:20/20;
K            = 1e-3;
ratios       = 0.1:0.2:0.9;
plies        = 20;
sigma_denoi  = 0;
% K                  = 0;
% for wavelength = [2 4 8 16 32 40]
process3 = process3.extract_local_orientation_3D_plywise_allwl(...
    PropertyName, wavelength, orientation, ratios, plies, K, sigma_denoi);

% show 3d slice
medfiltersize = [3 3 3];
angle_compens = -14; % degree
xslice        = 100 / process3.fx * 1e3;
yslice        = 100 / process3.fy * 1e3;
zslice        = 600 / process3.fs * 1e6; % 1.8 us
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

% fx_Scrollable_3d_view(abs(process3.img_hil));

%% save all figures
% "F:\Xiayang\results\Woven_samples\01032021_12h31m02s\f_5_omega07"
% "F:\Xiayang\results\Woven_samples\01032021_12h31m02s\f_10_omega07"
FolderName = "F:\Xiayang\results\Woven_samples\01032021_12h31m02s\f_10_omega07";   % the destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
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