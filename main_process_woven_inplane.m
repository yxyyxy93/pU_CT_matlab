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

load(strcat(FolderName, '20220428_BP1_3_025m_V313_25dbfullsize_PEmat.mat'));
% '20220207_BP1_3_025m_V313_25db_PEmat.mat'
% '20220303_BP1_3_025m_h5m_10db_PEmat.mat'
% '20220428_BP1_3_025m_V313_25db_filp_PEmat.mat'
% '20220428_BP1_3_025m_V313_25db_PEmat'
% '20220428_BP1_3_025m_V313_47db_PEmat'
% '20220428_BP1_3_025m_V313_25dbfullsize_PEmat'

%% EDA
% orthosliceViewer(process3.img);

process3.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500]; 
window   = [101 356 101 356 1 1500]; 
% window   = [171 270 1 100 1 1500];  % for 025m impact samples
% window   = [61 140 57 136 1 1500];  % for 025m impact samples

x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;

% % 3D viewer
orthosliceViewer(abs(process3.img_hil), "Colormap", jet);

%% show A scan
x   = 100;
y   = 100;

[Ascan, t_space, fss] = process3.demo_Ascan(x, y, strcat('A_scan_', Filename1(1:end-5)));

%% 2nd fft viewer
% p =temp = process3.img;
% tem fft(temp, [], 3);
% temp = abs(temp(:, :, 1:end/2));
% temp = fft(temp, [], 3);
% temp = abs(temp(:,:,1:end/2));
% fx_Scrollable_3d_view(temp);

%% autocorrelation
PropertyName = 'img_hil';
process3     = process3.caculate_autocorr(PropertyName);

orthosliceViewer(abs(process3.img_autoc), 'Colormap', jet);

%% surface calculatio - for autocorrelation
[lx, ly, ~]      = size(process3.img_autoc);
process3.front_I = ones(lx, ly);
process3.rear_I  = process3.front_I + 0.7e-6 * process3.fs; % 3.8 us;

%% normal time window

clc;
close all;

flag = 0;
delay       = 110;
max_len     = 600;
front_I_max = 500;

% PropertyName = 'img_hil_filter';
% PropertyName = 'img_hil_noise';
PropertyName = 'img_hil';
% PropertyName = 'img_WienerDeconv';
process3 = process3.find_damage_timewin(PropertyName, max_len, flag, delay, front_I_max, 1);

process3.show_surfaces('surface');


%% align the front surfaces
delay_I      = 100;
% PropertyName = 'img_hil';
% PropertyName = 'img';
process3     = process3.align_front_surface(delay_I, PropertyName);

process3.rear_I = process3.front_I + 200;
% back surface
% process3     = process3.align_rear_surface(delay_I, PropertyName);

%% show
% % close all;
% win_x = 1:100;
% win_y = 1:100;
% 
% process3.show_oneinterply_2dfft(7, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% 
% process3.show_oneinterply_2dfft(13, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply_2dfft(14, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply_2dfft(15, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% 
% close all;
% process3.show_oneinterply_2dfft(9, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply_2dfft(10, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply_2dfft(11, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% process3 = process3.show_oneinterply_2dfft(14, 'logGabor', 3000, win_x, win_y, TOF_oneply);

% %% save images 
% PropertyName = 'img_hil';
% ratios       = 0.2:0.1:0.8;
% sigma        = 5e-3;
% imagename    = 'C_scan_inam_plywise';
% 
% [process3, image1, ~]  = process3.define_plywise_inamCscan(1, ratios, PropertyName);
% process3.curvelet_denoise_inamCscan(sigma, imagename);
% [process3, image2, ~]  = process3.define_plywise_inamCscan(2, ratios, PropertyName);
% process3.curvelet_denoise_inamCscan(sigma, imagename);
% [process3, image3, ~]  = process3.define_plywise_inamCscan(3, ratios, PropertyName);
% process3.curvelet_denoise_inamCscan(sigma, imagename);
% [process3, image4, ~]  = process3.define_plywise_inamCscan(4, ratios, PropertyName);
% process3.curvelet_denoise_inamCscan(sigma, imagename);
% 
% % save('wovenstrcuture_5MHz.mat','image1','image2','image3','image4');
% 
% %% save images 
% imagename = 'C_scan_inam';
% 
% [process3, image1, ~]  = process3.define_parallel_inamCscan(1/20, PropertyName);
% process3.curvelet_denoise_inamCscan(sigma, imagename);
% [process3, image2, ~]  = process3.define_parallel_inamCscan(2/20, PropertyName);
% process3.curvelet_denoise_inamCscan(sigma, imagename);
% [process3, image3, ~]  = process3.define_parallel_inamCscan(3/20, PropertyName);
% process3.curvelet_denoise_inamCscan(sigma, imagename);
% [process3, image4, ~]  = process3.define_parallel_inamCscan(4/20, PropertyName);
% process3.curvelet_denoise_inamCscan(sigma, imagename);
% 
% save('wovenstrcuture_5MHz.mat','image1','image2','image3','image4');

%% read reference signal
% [Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
% filename               = strcat(Pathname1, Filename1);
% x                      = 10;
% y                      = 10;
% process3               = process3.read_refer(filename, x, y);
% 
% % cut the ref. signal
% win                    = [0.1e-6 2e-6]; % unit: s
% process3               = process3.cut_reference_signal(win);
% 
% % calculate the SNR
% win_noise              = [0.1e-6 0.5e-6];
% win_signal             = [0.6e-6 1e-6]; % unit: s
% snr                    = process3.calculate_SNR(win_signal, win_noise);
% 
% % align the signal
% process3               = process3.align_refer_ascan(x, y); 


%% one-plane EDA
% % 2D demo
% % define the C_scan_inam as well

%
% PropertyName = 'img_hil_filter';
% PropertyName = 'img_autoc';
% PropertyName = 'img_WienerDeconv';
% PropertyName = 'img_WienerDeconv_AR';
PropertyName = 'img_hil';
imagename    = 'C_scan_inam';

% ply    = 7;
% ratios = [0.2:0.1:0.8];
% ratios = 0.1;
% pename = 'C_scan_inam';
process3 = process3.define_parallel_inamCscan(1.5/20, PropertyName);
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
% % *************rocess3  = process3.define_plywise_inamCscan(ply, ratios, PropertyName);

% glabol RT
% image    = process3.C_scan_inam;
% [lx, ly] = size(image);
% r        = min(lx/2, ly/2) - 1;
% theta    = 1:180;
% [R, anguler_1D, ~, ~] = fx_Radonto1Dangular_correct(image, [lx/2+1, ly/2+1], r, theta);


% process3 = process3.show_oneinterply_2dfft(5, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% imagename   = 'mono_img';
% 2d log-Gabor fitler
wavelength  = (4:2:16) * 2;
orientation = 0:1:180; % degree
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

%% one-plane TV denoiser


%% one-plane EDA 2D structure tensor
% PropertyName = 'img_hil';
% 
% process3  = process3.define_parallel_inamCscan(4.5/24, PropertyName);
% ImageName = 'C_scan_inam';
% 
% sigma1 = 1;
% sigma2 = 1;
% process3.structural_tensor_2D(sigma1, sigma2, ImageName);

%% one-plane denoise 
close all;


process3 = process3.define_parallel_inamCscan(1.5/20, 'img_hil');
imagename   = 'C_scan_inam';
image1 = abs(process3.(imagename));

psd = abs(fftshift(fft2(image1))).^2;
figure, subplot(2,1,1);
imagesc(image1); 
% colormap jet; 
colorbar;
title('origianl image1');
subplot(2,1,2);
imagesc(psd); 
colormap jet; 
colorbar;
title('psd image1');

PSTDR = max(image1, [], 'all') / std(image1, 0,'all');
fprintf('PSTDR1 = %.2f \n', PSTDR);

%
outim=TVL1denoise(image1, 1, 100);
figure; 
imagesc(outim); 
colormap jet; 
colorbar;

%% one-plane EDA 2D analytic-signal
% PropertyName = 'img_autoc';
PropertyName = 'img_hil';
% PropertyName = 'img_WienerDeconv';
% PropertyName = 'img_WienerDeconv_AR'';

process3 = process3.define_parallel_inamCscan(4.5/20, PropertyName);
% process3.C_scan_inam = TVL1denoise(process3.C_scan_inam , 1, 100);
imagename   = 'C_scan_inam';
clc;
close all;

s_0    = 5;
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

mask_size = s_0*2;
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
% 
% orien_map(apexa_map>=0) = nan;

%% 'distance to front and rear' in-plane orientation extraction 2D analytic-signal 
% 2D analytic-signal
% PropertyName  = 'img_hil';
PropertyName  = 'img_autoc';
% PropertyName  = 'img_hil_filter';
% PropertyName = 'img_WienerDeconv';

ratios        = 0/20:0.1/20:20/20;
sigma_denoise = 0;
% 2D analytic signal

tic;
process3 = process3.extract_local_orientation_3D_parallel_2Danalytic_signal...
    (PropertyName, s_f, s_c, mask_size, ratios, sigma_denoise);    
toc;

%% show slice
angle_compens = 14;
% volumeViewer(process3.Inplane_direction_3D_ID + angle_compens);

process3.statistic_angular_distribution(angle_compens);
caxis([0 0.02]);

%%  glabol RT - 3D
ratios = 0/24:0.1/24:24/24;
ImageName = 'C_scan_inam';
image    = process3.C_scan_inam;
[lx, ly] = size(image);
r        = min(lx/2, ly/2) - 1;
theta    = 1:180;
anguler_3D = zeros(length(ratios), length(theta));
for i = 1:length(ratios)
    PropertyName = 'img_hil';
    process3  = process3.define_parallel_inamCscan(ratios(i), PropertyName);
    [R, anguler_1D, ~, ~] = fx_Radonto1Dangular_correct(process3.C_scan_inam,...
        [lx/2+1, ly/2+1], r, theta);
    anguler_3D(i, :) = anguler_1D/sum(anguler_1D);
    clc;
    fprintf('RT progress: %0.2f%%\n',100*i/length(ratios));
end

%% monogenic signal analysis
% PropertyName = 'img_hil';
% process3     = process3.define_parallel_inamCscan(3.5/20, PropertyName);
% imagename    = 'C_scan_inam';
% % fillna
% C_scan_inam = process3.(imagename);
% C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
% 
% % ******** get the image to negative and positive periodic
% C_scan_inam = C_scan_inam - mean(C_scan_inam);
% 
% %
% [Y, X]       = size(C_scan_inam);
% cw           = 10;
% filtStruct   = createMonogenicFilters(Y, X, cw, 'lg', 0.55);
% [m1, m2, m3] = monogenicSignal(C_scan_inam, filtStruct);
% % [m1, m2, m3] = monogenicSignal_nofilter(C_scan_inam, filtStruct);
% 
% % close all;
% % Local energy (calculated on a per-scale basis)
% LE = localEnergy(m1, m2, m3);
% % Local phase (calculated on a per-scale basis)
% LP = localPhase(m1, m2, m3);
% % Local orientation (calculated on a per-scale basis)
% % Only need to pass the odd parts (m2,m3) as even part (m1) is irrelevant
% LO = localOrientation(m2, m3);
% %
% close all;
% % Display
% figure;
% imagesc(C_scan_inam), axis image; colormap(jet);
% h = colorbar;
% set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
% % colormap gray
% title('Test Image');
% 
% figure();
% imagesc(reshape(LE, Y, [])), axis image;
% h = colorbar;
% colormap(jet);
% set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
% title('Local Energy Over Scales');
% figure();
% imagesc(reshape(LP, Y, [])), axis image, axis off, colormap gray;
% title('Local Phase Over Scales');
% figure();
% imagesc(reshape(LO, Y, [])), axis image, axis off; colorbar;
% title('Local Orientation Over Scales (radians)');
% h = colorbar;
% colormap(hsv);

%% 3D monogenic signal analysis
% clc;
% D         = process3.img;
% [Y, X, Z] = size(D);
% 
% % woven strucutre
% cw_xy = 10;
% cw_z  = 250e6/12.8e6;
% % D     = imgaussfilt3(D, 1);
% % Now use these wavelengths to create a structure containing
% % frequency-domain filters to calculate the monogenic signal. Same as for
% % the 2D case, except now we need to pass three image dimensions and the
% % two sets of wavelengths
% filtStruct = createMonogenicFilters3D(Y, X, Z, cw_xy, cw_z, 'lg', 0.55);
% % Now we can use this structure to find the monogenic signal for the volume
% [m1,m2,m3,m4] = monogenicSignal3D(D, filtStruct);
% % Local energy (calculated on a per-scale basis)
% LE = localEnergy3D(m1,m2,m3,m4);
% % Local amplitude (calculated on a per-scale basis)
% LA = sqrt(LE);
% % Local phase (calculated on a per-scale basis)
% LP = localPhase3D(m1,m2,m3,m4);
% [FS,FA] = featureSymmetry3D(m1,m2,m3,m4);
% 
% % PC = phaseCongruency3D(m1,m2,m3,m4,0.05);
% 
% % % Cscan image
% % z  = 400;
% % im = LA(:, :, z);
% % imagesc(im), axis image; colormap(jet);
% % h = colorbar;
% % set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
% % 
% % % Ascan
% % x = 45;
% % y = 45;
% % figure,plot(squeeze(LA(x, y, :)));
% 
% LP_1D = abs(process3.img_hil);
% 
% fx_Scrollable_3d_view(LP_1D);

%% multi one-plane EDA 
% % 2D demo
% ply_num             = 20;
% ply_wise_Cscan_inam = cell(1, ply_num);
% Idof_N_arr          = NaN(ply_num, 180);
% % PropertyName        = 'img_hil';
% % PropertyName        = 'img_hil_filter';
% PropertyName        = 'img_WienerDeconv';
% % imagename           = 'C_scan_inam_plywise';
% imagename           = 'C_scan_inam';
% 
% for ply = 1:ply_num
% %     process3     = process3.show_oneinterply_2dfft(ply, 'logGabor', 3000, win_x, win_y, TOF_oneply);
% %     process3 = process3.define_plywise_inamCscan(ply, ratios, PropertyName);
%     process3 = process3.define_parallel_inamCscan(ply/2/ply_num, PropertyName);
%     process3 = process3.compute_logGabor_filter_withoutFig( ...
%         PropertyName, wavelength, orientation, SFB, SAR, imagename);
%     [image_ori, Idof_N] = process3.show_orientation_by_ID_allwl( ...
%         wavelength, orientation, K, imagename, angle_compens);
%     ply_wise_Cscan_inam{ply} = image_ori;
%     Idof_N_arr(ply, :) = Idof_N;
%     close all;
% end
% 
% % display all images in on figure
% close all;
% figure('Name', 'ply_wise_Cscan_inam')
% for ply = 1:ply_num
%     subplot(ply_num/10, 10, ply)
%     pcolor(ply_wise_Cscan_inam{ply}); shading flat;
%     colormap hsv;
%     caxis([-90 90]);
%     title(num2str(ply));
% end
% colorbar;
% 
% % **** 
% Stacking_sequence = [
%     45 0 -45 -90 ...
%     45 0 -45 -90 ...
%     45 0 0 45 ...
%     -90 -45 0 45 ...
%     -90 -45 0 45];
% maps = cell(1, ply_num);
% for i = 1: ply_num
%     maps{i} = ['Ply' num2str(i) ':    ' num2str(Stacking_sequence(i)) '\circ'];
% end
% %
% Idof_N_arr = circshift(Idof_N_arr, 22, 2);
% fx_creat_chartPlot(maps, Idof_N_arr, Stacking_sequence);
% %*********
% 
% %
% figure('Name', 'ply_wise_anguler_1D')
% anguler_stacking = nan(180, ply_num);
% for ply = 1:ply_num
%     subplot(ply_num/4, 4, ply)
%     plot(-89:90, anguler_1D_ply_wise{ply}, 'linewidth', 2);
%     anguler_stacking(:, ply) = anguler_1D_ply_wise{ply};
% end
% 
% fx_Scrollable_3d_view(abs(process3.img_hil));

%% 'distance to front and rear' in-plane orientation extraction
% % RT
% PropertyName  = 'img_hil';
% % PropertyName  = 'img_hil_filter';
% % PropertyName = 'img_WienerDeconv';
% 
% % radius        = 10;
% % theta         = 1:1:180; % degree
% ratio         = 0/20:0.1/20:20/20;
% % sigma_denoise = 0;
% % tic;
% % process3      = process3.extract_local_orientation_RT_3D_parallel(...
% %     PropertyName, radius, theta, ratio, sigma_denoise);
% % toc;
% % 
% % % show slice
% % xslice        = 100 / process3.fx * 1e3;
% % yslice        = 100 / process3.fy * 1e3;
% % zslice        = []; % us
% % angle_compens = 0;
% % process3.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);
% 
% % 2D log-Gabor filter
% wavelength  = 10:2:10;
% orientation = 0:1:180; % degree
% SFB         = 1.5;       % [0.5 2.5]
% SAR         = 0.5;     % [0.23 0.92]
% K           = 0;   
% % sigma              = 5e-4;
% sigma       = 0;
% tic;
% process3 = process3.extract_local_orientation_3D_parallel_allwl(...
%     PropertyName, wavelength, orientation, ratio, K, sigma);
% toc;
% 
% 
% % show slice
% xslice        = 50 / process3.fx * 1e3;
% yslice        = 50 / process3.fy * 1e3;
% zslice        = [];
% mfsize        = [1 1 1];
% angle_compens = -14;
% process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

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
    (PropertyName, wavelength, orientation, z_theta, z_range, K, SFB, SAR);
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