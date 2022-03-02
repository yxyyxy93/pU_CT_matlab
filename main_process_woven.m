% read the file
% used to estimate thickness as well
% This is used for the test the orignal phase-derived interply track and the log-gabor fitelr
% for 1st paper parametric study

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
FolderName = "F:\Xiayang\results\Woven_samples\";   % the destination folder
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '20220207_BP1_3_025m_V313_25db_PEmat.mat'));
% '01032021_12h31m02s.mat'
% '01032021_14h17m04s.mat'
% '01032021_15h41m01s.mat'

%% EDA
process3.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500]; 
% window = [1 120 151 300 1 1500]; 
window   = [1 400 1 400 200 1200];  
x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;
% 3D viewer
volumeViewer(abs(process3.img_hil));

%% show A scan
% define the index to select Ascan 
x           = 45;
y           = 45;
process3.demo_AS_3D_inclinfq(x, y);

%% show C scan
% define the index to select Ascan 
z            = 400;
PropertyName = 'img_hil';
process3.show_Cscan(z, PropertyName);

%%

process3.img_hil_filter = LA;
PropertyName = 'img_hil_filter';
% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y

MinPD   = 25;
MinPH   = 0.015;  % these 2 parameters need to be changed for surface estimation.
% surface calculation
max_len = 1300;
alpha   = 4e-3;
A_ratio = 0.5;

process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, 40, 150);

% A_ratio      = 0.4; % for 0.3 m drop
process3     = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
process3.show_surfaces(filename_fig(1:end-5));


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
    85, 146; % x
    50, 300; % y
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
disp(' ')
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
index = 150;
Bwin  = 1:400;
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
zslice = 250; % near the middle
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
xslice = 150; % near the middle
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

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
MinPD        = 25;
MinPH        = 0.03;  % these 2 parameters need to be changed for surface estimation.
PropertyName = 'img_hil';
x = 22;
y = 22;
process3.show_Ascan_inam_peaks(x, y, MinPD, MinPH, PropertyName); % x, y

alpha   = 4e-3;
A_ratio = 0.8;
max_len = 1000;
x       = 30;
y       = 290;
% A_ratio      = 0.4; % for 0.3 m drop
process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, x, y)
        
process3     = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
process3.show_surfaces(filename_fig(1:end-5));

process3     = process3.recover_surface;
process3     = process3.smooth_rear_I;
process3.show_surfaces(filename_fig(1:end-5));

TOF_walls    = -mean(process3.front_I-process3.rear_I, 'all', 'omitnan')/process3.fs;
disp(TOF_walls);
% %
% process      = process.smooth_rear_I;
% process.show_surfaces(filename_fig(1:end-4));
% %
% process      = process.recover_surface;

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

%% save all figures"z
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