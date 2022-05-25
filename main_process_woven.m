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

% % read the settings from excel
% [Filename1, Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
% settinig_file          = strcat(Pathname1, Filename1);
% 
% % read the data from the struct, or reset the dataset
% process3               = process3.loadSettings(settinig_file);

process3               = process3.read_origin_data; % read (reset) the dataset

%% ***********
FolderName = "C:\Users\xiayang\OneDrive - UGent\matlab_work\results\Woven_samples\";   % the destination folder
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '20220303_BP1_3_025m_h5m_10db_PEmat.mat'));
% '20220207_BP1_3_025m_V313_25db_PEmat.mat'
% '01032021_14h17m04s.mat'
% '01032021_15h41m01s.mat'

%% EDA
% orthosliceViewer(process3.img);

process3.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500]; 
% window   = [200 350 270 420 300 1300]; 
window   = [1 270 321 500 1 1500];  
x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;
% % 3D viewer
% volumeViewer(abs(process3.img_hil));

%% show A scan
% define the index to select Ascan 
x   = 35;
y   = 35;
% process3.demo_AS_3D_inclinfq(x, y);

process3.show_hilbert_Ascan(x, y);

%% show C scan
% define the index to select Ascan 
close all;
for z = 1000:50:1100
    % z = ;
    PropertyName = 'img';
    process3 = process3.show_Cscan(z, PropertyName);
%     process3.check2dfft_inclphase(z, PropertyName); % 2d spectrum
end

%% low-pass filter
bandpassFreq = 30e6;
process3     = process3.Filter_lowpass(bandpassFreq, 'img_hil');

%% surface search
global min_pks 
min_pks      = 0.5;
PropertyName = 'img_hil_filter';
% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y
clc;
% close all;

MinPD   = 30;
MinPH   = 0.015;  % these 2 parameters need to be changed for surface estimation.
% surface calculation
max_len = 1499;
alpha   = 4.5e-3;
A_ratio = 0.99;

process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, 220, 55);

%%
% A_ratio      = 0.4; % for 0.3 m drop
process3     = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
close all;
process3.show_surfaces(filename_fig(1:end-5));

%% check 2 times fft
temp = process3.img;
temp = temp(209, 277, :);
temp = squeeze(temp);

figure, subplot(3, 1, 1);
plot(temp);
subplot(3, 1, 2);
temp_fft = fft(temp, 2^14);
plot(abs(temp_fft));
subplot(3, 1, 3);
temp_fft2 = fft(abs(temp_fft));
plot(abs(temp_fft2(1:end/2)));

%%
PropertyName = 'img_hil';
max_len      = 1500;
flag_DG      = 1; % apply DynamicGate
delay        = 50; % delay for the search on the 2nd fft, unit: points
filter_flag  = 1;
passband     = [0.5e6 15e6];

process3 = process3.find_surface_2fft(PropertyName, max_len, flag_DG, delay, filter_flag, passband);

close all;
process3.show_surfaces(filename_fig(1:end-5));

%% internal damage features
process3.damage_imaging;

%% fingerprint one ply - sino surface fitting
% z        = 200;
% z_range  = 100:10:1000;
z_range  = 5/20:5/20:5/20;
best_object_3D   = nan(length(z_range), 1); 
best_location_3D = nan(length(z_range), 2); 
for z_index = 1:length(z_range)
    z = z_range(z_index);
    [best_object_z, best_location_z] = process3.fingerprint_oneply(z, 'img_hil');
    best_object_3D(z_index)      = best_object_z;
    best_location_3D(z_index, :) = best_location_z;
    clc;
    disp(['Procedure: ', num2str(100*z_index/length(z_range)), '%']);
end

%% make fingerprint along z depth
ratios       = 0/20:0.5/20:20/20;
PropertyName = 'img_hil';
% process3     = process3.fingerprint_2dfft_fitsurface(ratios, PropertyName);
process3     = process3.fingerprint_2dfft_fitsurface_v2(ratios, PropertyName);

%%
temp  = process3.fft2d_pha_fingerprint;
temp1 = squeeze(sum(temp, 1, 'omitnan'));
temp1 = temp1.';
figure,
imagesc(temp1);

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

%%
% out-of-plane angles
f0    = 6.3e6;
sigma = 0.6;
process3 = process3.Filter_logGabor(f0, sigma, 'img');

% smoothing scale
ds_rate      = 1;
sigma1       = [3, 3, 3];
% integration scale
sigma2       = [3, 3, 3];
PropertyName = 'img_hil_filter';
process3     = process3.structural_tensor(sigma1, sigma2, PropertyName, ds_rate);

% display the angles 3d slices
medf_kernel  = [3, 3, 3];
xslice       = 100 / process3.fx * 1e3;
yslice       = 100 / process3.fy * 1e3;
zslice       = [];
process3.show_angles_ST_zangle(medf_kernel, xslice, yslice, zslice, ds_rate);

%
B_type     = 'x';
index      = 100;
Bwin       = 1:150;
process3.show_angles_ST_Bscan(B_type, index, Bwin);

%
B_type     = 'y';
index      = 100;
Bwin       = 1:150;
process3.show_angles_ST_Bscan(B_type, index, Bwin);

temp = real(process3.angle_x);
temp(isnan(temp)) = 0;
temp              = temp / pi * 180;
orthosliceViewer(temp, 'ScaleFactors', [1 1 0.2], 'colormap', jet, 'DisplayRange', [-2 2]);

%% 3d ply track

% % log-Gabor filter
% f0        = 6.3e6;
% sigma     = 0.7;
% process3  = process3.Filter_logGabor(f0, sigma, 'img_hil');

% low-pass filter
bandpassFreq = 10e6;
process3     = process3.Filter_lowpass(bandpassFreq, 'img_hil');

threshold    = 0.01;
process3     = process3.track_interply('img_hil_filter');
% process3  = process3.track_interply_inph(threshold, 'img_hil_filter', 21);

% use 2nd-harmonic for first and the last interplies, 
% use fundamental resonance for another interplies
f0_1     = 12.8e6;
sigma0_1 = 0.7;
f0_2     = 6.3e6;
sigma0_2 = 0.6;
nol      = 21; % make it 25+1 in case that there is redundant accidently-tracked interply
process3 = process3.track_interply_hybrid( ...
    'img', f0_1, sigma0_1, f0_2, sigma0_2, nol);

%%
TOF_walls    = -mean(process3.front_I - process3.rear_I, 'all', 'omitnan');

close all;
B_type     = 'x';
index      = 80;
Bwin       = 1:270;
TOF_oneply = TOF_walls / 20; % 24 plies
process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

% process.show_track_interply(xslice, yslice, zslice);
win_x = 1:270;
win_y = 1:100;
process3.show_oneinterply(3, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(10, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(15, 'logGabor', 3000, win_x, win_y, TOF_oneply);

%%
close all;
win_x = 1:270;
win_y = 1:100;
process3.show_oneinterply_2dfft(1, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply_2dfft(2, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply_2dfft(3, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply_2dfft(4, 'logGabor', 3000, win_x, win_y, TOF_oneply);

%% knn
process3.show_interply_track_3D;


%% apply 2d fft filter
close all;

PropertyName = 'img_hil';
% filtertype   = {'BP', 'Butterworth'};
filtertype   = {'norchBP', 'Butterworth'};
slicetype    = 'surface parallel';
z            = 5/20;
process3.apply2dfft_filter(z, PropertyName, filtertype, slicetype)

%% amplitude rise method to determine the woven structure - spatial filter
clc;
close all;
PropertyName = 'img_hil';
zrange       = [450 500];
drop         = 0.1; % dB
sigma        = 1;

process3 = process3.amplitude_rise_oneimage_spatialfilter(PropertyName, zrange, drop, sigma);

%% 2d slice and 2d fft
% 
z            = 450;
PropertyName = 'img_hil';
% 
% process3.show_Cscan(z, PropertyName);

% process3.show_Cscan_2dfft(300, PropertyName);

% roi = [
%     75, 200; % x
%     100, 300; % y
%     ]; % for 20220207_BP1_3_025m_V313_25db_PE

roi = [
    1, 150; % x
    1, 150; % y
    ]; % for 20220207_BP1_3_025m_V320_25db_PE

% 2d fft filter
% LP
close all;
R_lo_H = 0;
R_Lo   = 15; 
process3.apply2dfft_LPandHPfilter(z, PropertyName, R_lo_H, R_Lo, roi);

% HP
close all;
R_lo_H = 60;
R_Lo   = 15; 
process3. apply2dfft_HPfilter(z, PropertyName, R_lo_H, roi);

%% 2d sinc function

% 2d fft filter
% 2d sinc funciotn
close all;

freq   = 1/30;
process3.apply2dfft_LPsinc(z, PropertyName, freq, roi);

%% 
% bandreject
close all;
D_0    = 50;
W      = 15;
n      = 2;
F_type = 'ideal';
process3.apply2dfft_bandreject(z, PropertyName, D_0, W, n, F_type, roi);
F_type = 'Butterworth';
process3.apply2dfft_bandreject(z, PropertyName, D_0, W, n, F_type, roi);
F_type = 'Gaussian';
process3.apply2dfft_bandreject(z, PropertyName, D_0, W, n, F_type, roi);

%% 
% band pass
close all;
D_0    = 25;
W      = 10;
n      = 2;
F_type = 'ideal';
process3.apply2dfft_bandpass(z, PropertyName, D_0, W, n, F_type, roi);
F_type = 'Butterworth';
process3.apply2dfft_bandpass(z, PropertyName, D_0, W, n, F_type, roi);
F_type = 'Gaussian';
process3.apply2dfft_bandpass(z, PropertyName, D_0, W, n, F_type, roi);


%%
% notch filter
close all;
D_0    = 10;
X_0    = [50 0 35 35];
Y_0    = [0 50 35 -35];
% X_0    = x/10;
% Y_0    = y/10;
n      = 2;
F_type = 'ideal';
process3.apply2dfft_notch(z, PropertyName, D_0, X_0, Y_0, n, F_type, roi);
D_0    = 30;
F_type = 'Butterworth';
process3.apply2dfft_notch(z, PropertyName, D_0, X_0, Y_0, n, F_type, roi);
D_0    = 30;
F_type = 'Gaussian';
process3.apply2dfft_notch(z, PropertyName, D_0, X_0, Y_0, n, F_type, roi);

%% A-contrario quasi-periodic noise removal
z            = 450;
PropertyName = 'img_hil';
% fillna
im_gt = abs(process3.(PropertyName)(:, :, z));
im_gt(isnan(im_gt)) = mean(im_gt, 'all', 'omitnan');

close all;
im_noise_remov = ACARPENOS(im_gt, 128, 1e0, 0, roi);

%% create multiple pathces
% calculate minimum power spectrum
z            = 450;
PropertyName = 'img_hil';
% fillna
im_gt = abs(process3.(PropertyName)(:, :, z));
%
size_patches = 128;
step_patches =round(size_patches/2); % translation between patches (in pixels)
[rows, cols] = size(im);
P = numel(1:step_patches:(rows-size_patches))*numel(1:step_patches:(cols-size_patches));  % total number of patches
disp(' ');
%
H           = hanning(size_patches)*hanning(size_patches)';  % patches are weighted by a Hann window to remove border effects
% H = 
logspeclist = zeros(size_patches,size_patches,P);
imagelsit   = zeros(size_patches,size_patches,P);

n = 0;
for i=1:step_patches:(rows-size_patches)
  for j=1:step_patches:(cols-size_patches)
  n                  = n+1;
  imagette           = H.*im_gt(i:i+size_patches-1,j:j+size_patches-1);
%   imagette           = im_gt(i:i+size_patches-1,j:j+size_patches-1);
  imagelsit(:,:,n)   = imagette;
  logspeclist(:,:,n) = 2*log(abs(fftshift(fft2(imagette))));
  end
end

% close all;
clc;
figure,
prow = 5;
pcol = 3;
for i = 1: prow*pcol
    subplot(prow, pcol, i);
    imagesc(imagelsit(:,:,i));
    axis image;
    colormap(jet);
end   

figure,
for i = 1: prow*pcol
    subplot(prow, pcol, i);
    imagesc(logspeclist(:,:,i));
    axis image;
    colormap(jet);
end  

%% 2D curvelet
z            = 450;
PropertyName = 'img_hil';
% fillna
noisy_img    = abs(process3.(PropertyName)(:, :, z));

disp('Take curvelet transform: fdct_usfft');
noisy_img(isnan(noisy_img)) = mean(noisy_img, [1 2], 'omitnan');
tic; C     = fdct_wrapping(double(noisy_img), 0, 1); toc;
% upsample the Coefficients by duplication
% save as Information diagram
[img_m, img_n] = size(noisy_img);
ori_len_max    = length(C{end-1});
C_temp         = NaN(img_m, img_n, length(C), ori_len_max);
for i = 1:length(C)
    for j = 1:length(C{i})
        curvelet_C        = abs(C{i}{j});  % absolute value in case of complex
        % if the C is empty in last scale
        if isempty(curvelet_C)
            break;
        end
        len_ori           = ori_len_max / length(C{i});
        % using bulid-in 'imresize' function
        curvelet_C_resize = imresize(curvelet_C, [img_m img_n]);
        %                     figure('Name', [num2str(i), '_' num2str(j)]), imagesc(curvelet_C_resize);
        for idx_len_ori = 1:len_ori
            C_temp(:, :, i, (j-1)*len_ori+idx_len_ori) ...
                = curvelet_C_resize / 2^(-3*i/2); % theoretical scale
            % In other words, one would need to compute on the order of 2^2j coefficients
            % per scale and angle as opposed to only about 2^3j/2 in the USFFT-based implementation.
        end
    end
end

clc;
figure,
prow = 8;
pcol = 8;
% plot
i = 6;
for j = 1:length(C{i})
    curvelet_C        = abs(C{i}{j});  % absolute value in case of complex
    subplot(prow, pcol, j);
    imagesc(curvelet_C); colormap(jet)
end

%% curvelet denoise
z            = 450;
PropertyName = 'img_hil';
% fillna
sigma        = 1e-1;

noisy_img    = abs(process3.(PropertyName)(:, :, z));
%             closeNaN_part               = imopen(NaN_part, se); % Morphologically close imag
restored_img           = fx_curvelet_denoise_enhanced(sigma, noisy_img);
close all
%             title('Original Image')
cf      = figure('Name', 'x-y slice-plane of instantaneous amplitude');
set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');

%
xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
%             figure; clf; imagesc(restored_img); colormap gray; axis('image');
%             title('Restored image');
C_scan_inam_denoise = restored_img;

% plot 2d fft filtered
ax     = subplot(1, 1, 1);
imagesc(abs(C_scan_inam_denoise));
axis image; colormap(jet);
hold on;
h = colorbar;
set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
set(ax, 'linewidth', 2);
% display snr
disp([F_type 'snr and cnr image_filter:']);
[snr, cnr] = fx_image_snr_2(C_scan_inam_denoise, roi);
disp([snr cnr]);

%% monogenic signal analysis
z            = 450;
PropertyName = 'img_hil';
% fillna
C_scan_inam = abs(process3.(PropertyName)(:, :, z));
C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');

% ******** get the image to negative and positive periodic
C_scan_inam = C_scan_inam - mean(C_scan_inam);
Ascan      = C_scan_inam(:,100,:);
Ascan_inam = abs(hilbert(Ascan));
plot(Ascan, 'linewidth', 2, 'displayname', 'signal');
hold on;
plot(Ascan_inam, 'linewidth', 2, 'displayname', 'instantaneous amplitude');
xlabel('Y displacement');
ylabel('Amplitude (arb.)')
set(gca, 'Fontname', 'times new Roman', 'Fontsize', 12);
set(gca, 'linewidth', 2) ;
legend;

% 
[Y, X]       = size(C_scan_inam);
cw           = [5 10 20 40];
filtStruct   = createMonogenicFilters(Y, X, cw, 'lg', 0.66);
[m1, m2, m3] = monogenicSignal(C_scan_inam, filtStruct);
% [m1, m2, m3] = monogenicSignal_nofilter(C_scan_inam, filtStruct);

close all;
% Local energy (calculated on a per-scale basis)
LE = localEnergy(m1, m2, m3);
% Local phase (calculated on a per-scale basis)
LP = localPhase(m1, m2, m3);
% Local orientation (calculated on a per-scale basis)
% Only need to pass the odd parts (m2,m3) as even part (m1) is irrelevant
LO = localOrientation(m2, m3);
% Display
cf = figure();
imagesc(C_scan_inam), axis image; colormap(jet);
h = colorbar;
set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
% colormap gray
title('Test Image');
%
close all
for i = 1:size(LE, 4)
    figure();
    C_scan_inam_denoise = sqrt(LE(:, :, :, i)); % local amplitdue = sqrt(local energy)! 
    imagesc(C_scan_inam_denoise), axis image;
    colormap jet;
    h = colorbar;
    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
    disp(['local energy, snr and cnr mono_' num2str(cw(i)) ':']);
    [snr, cnr] = fx_image_snr_2(C_scan_inam_denoise, roi);
    disp([snr cnr]);
end

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

%% ***************** 3D analysis below *************************** %
% Bscan origin
B_type = 'y';
index = 100;
Bwin  = 1:150;
PropertyName = 'img_hil';
process3.demo_Bscan_inst(B_type, index, Bwin, PropertyName);

%% 3D monogenic signal analysis
clc;
D         = process3.img;
[Y, X, Z] = size(D);

cw_xy = 5;
cw_z  = 250e6/15e6;

% woven strucutre
cw_xy = 30;
cw_z  = 250e6/15e6;

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

% Cscan image
z  = 400;
im = LA(:, :, z);
imagesc(im), axis image; colormap(jet);
h = colorbar;
set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');

% Ascan
x = 45;
y = 45;
figure,plot(squeeze(LA(x, y, :)));

volumeViewer(LA);

%% 3D monogenic signal display
% Display one slice
% 
zslice = 50; % near the middle
figure()
imshow(D(:,:,zslice)), axis image, colormap jet
title('Test Volume Slice')

figure()
imagesc(reshape(LE(:,:,zslice,:),Y,[])), axis image, colormap jet
title('Local Energy Over Scales')

figure()
imagesc(reshape(LP(:,:,zslice,:),Y,[])), axis image, colormap gray
title('Local Phase Over Scales')

figure()
imagesc([FS(:,:,zslice),FA(:,:,zslice)]), axis image, colormap jet
title('3D Feature Symmetry and Asymmetry');

% 
close all;
xslice = 50; % near the middle
figure()
imgshow = squeeze(D(:,xslice,:));
imshow(rot90(imgshow, 3)), axis image, colormap jet
title('Test Volume Slice')

figure()
imgshow = reshape(LE(:,xslice,:,:),Y,[]);
imagesc(rot90(imgshow, 3)), axis image, colormap jet
title('Local Energy Over Scales')

figure()
imgshow = reshape(LP(:,xslice,:,:),Y,[]);
imagesc(rot90(imgshow, 3)), axis image, colormap gray
title('Local Phase Over Scales')

figure()
imgshow = [squeeze(FS(:,xslice,:)),squeeze(FA(:,xslice,:))];
imagesc(rot90(imgshow, 3)), axis image, colormap jet
title('3D Feature Symmetry and Asymmetry');

%% 3D fft filter 

process3            = process3.nDfft;
[lxwn, lywn, lfreq] = size(process3.img_3dfft);

%
sum_fx = squeeze(sum(abs(process3.img_3dfft, 2)));
sum_yf = squeeze(sum(abs(process3.img_3dfft, 2)));
sum_xy = squeeze(sum(abs(process3.img_3dfft, 3)));

figure,
[x, y] = meshgrid(1:lxwn, 1:lywn);
surf(x, y, ones(lxwn, lywn), abs(sum_xy));
hold on;

surf(x, 1, z, abs(sum_fx));
hold on;
surf(1, y, z, abs(sum_yf));
hold on;

%%
volumeViewer(abs(process3.img_3dfft));
sliceViewer(abs(process3.img_3dfft));

%% original results display
xslice = 50 / process3.fx * 1e3;
yslice = 50 / process3.fy * 1e3;
zslice = []; % us

process3.show_inaminph_3D(xslice, yslice, zslice);

z            = 400;
PropertyName = 'img_hil';
process3.show_Cscan(z, PropertyName);

%% log-Gabor filter
% % filtered
f0        = 5e6;
sigma     = 0.7;

% for debug
% threshold = 0.05;
x         = 40;
y         = 40;
% process3.demo_logGabor_plytrack_inph_v2(x, y, f0, sigma, threshold);
process3.show_hilbert_Ascan(x, y);
% process3.show_logGabor_Scaleogram(3e6:0.1e6:16e6, sigma, x, y);

% 0.676 us - 2.13 us
%  
process3 = process3.Filter_logGabor(f0, sigma, 'img');
nol      = 49;

%% Volume analysis
ratio        = 0.5;
PropertyName = 'img_hil';
[process3, C_scan_inam_temp, C_scan_index_temp] = process3.define_parallel_inamCscan(ratio, PropertyName);
imagesc(C_scan_index_temp);

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