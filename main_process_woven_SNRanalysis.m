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
% '20220428_BP1_3_025m_V313_25dbfullsize_PEmat'

%% EDA
% orthosliceViewer(process3.img);

process3.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500]; 
window   = [1 128 1 128 1 1500]; 
% window   = [171 270 1 100 1 1500];  % for 025m impact samples
% window   = [61 140 57 136 1 1500];  % for 025m impact samples

x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;

% % 3D viewer
orthosliceViewer(abs(process3.img_hil), "Colormap", jet);

%% show A scan
close all;

x   = 100;
y   = 100;

[Ascan, t_space, fss] = process3.demo_Ascan(x, y, strcat('A_scan_', Filename1(1:end-5)));

process3.show_hilbert_Ascan(x, y);

figure;
stft(real(Ascan), process3.fs);

%% normal time window
PropertyName = 'img_hil';
% PropertyName = 'img_hil_filter';
delay        = 600;
max_len      = 1200;
flag         = 0;
front_I_max  = 570;
process3     = process3.find_damage_timewin_asignsurface(PropertyName, max_len, flag, delay, front_I_max, 1);
process3.show_surfaces(filename(1:end-5));

delay_I      = 300;
% PropertyName = 'img_hil';
% PropertyName = 'img';
% process3     = process3.align_front_surface(delay_I, PropertyName);

%% extract ascan and calcualte snr separately

x = 20;
y = 20;

Ascan = real(squeeze(process3.img_hil(x, y, :)));
front_point = process3.front_I(x, y);
rear_point = process3.rear_I(x, y);

figure,
plot(Ascan);
hold on;
total_point = 21;
points_arr = round(linspace(front_point, rear_point, total_point));
scatter(points_arr, Ascan(points_arr), 5, 'red','filled');

%
img_3D = abs(process3.img_hil);
% orthosliceViewer(img_3D);

Ascan_noise = img_3D(:, :, 1:200);
noise_power = mean(Ascan_noise.^2, "all");

front_ave = mean(process3.front_I, "all");
rear_ave  = mean(process3.rear_I, "all");
points_ave_arr = round(linspace(front_ave, rear_ave, total_point));

% SNR calculation
snr_arr = nan(total_point-1, 1);
for i = 1:total_point-1
    Ascan_sub    = img_3D(:,:,points_ave_arr(i):points_ave_arr(i+1));
    signal_power = mean(Ascan_sub.^2, "all");
    snr_arr(i)   = 10*log10(signal_power/noise_power);
end

figure,
plot(snr_arr); 
xlabel('depth range (from front surface to back surface)');
ylabel('SNR (dB)');

%% read reference signal
[Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
filename               = strcat(Pathname1, Filename1);
x        = 10;
y        = 10;
process3 = process3.read_refer(filename, x, y);

% cut the ref. signal
win      = [0.3e-6 1.9e-6]; % unit: s
process3 = process3.cut_reference_signal(win);

% align the signal
process3 = process3.align_refer_ascan(x, y); 

% process3.show_reference_signal;
%% show A scan - different center frequencies
x   = 10;
y   = 75;
Fcs = [10, 5] * 1e6;
snr = 30;
process3.show_Ascan_freqManipulation(x, y, Fcs, snr);

%% frequency manipulation dataset
clc;
close all;
Fc            = 5e6;
process3      = process3.data_freqManipulation(Fc);
sampling_rate = 1; % for reset the sampling_rate

img_3D_fc = abs(process3.img_hil_noise);

% plot
figure,
subplot(2,1,1);
image1 = img_3D(:,:,points_ave_arr(4));
imagesc(image1); colormap jet; colorbar;
subplot(2,1,2);
image2 = img_3D_fc(:,:,points_ave_arr(4));
imagesc(image2); colormap jet; colorbar;


figure,
subplot(2,1,1);
imagesc(img_3D(:,:,points_ave_arr(15))); colormap jet; colorbar;
subplot(2,1,2);
imagesc(img_3D_fc(:,:,points_ave_arr(15))); colormap jet; colorbar;


%%
% add noise to lower the snr
[lx, ly, lz]   = size(img_3D);
snr            = 25;
sigPower       = sum(abs(img_3D).^2, 3) / lz; % linear
reqSNR         = 10^(snr/10);
noisePower     = sigPower/reqSNR;
noise          = repmat(sqrt(noisePower), [1 1 lz]) .* randn(lx, ly, lz);
img_noiseadded = img_3D + noise;


Ascan_noise = img_noiseadded(:, :, 1:200);
noise_power = mean(Ascan_noise.^2, "all");

front_ave = mean(process3.front_I, "all");
rear_ave  = mean(process3.rear_I, "all");
points_ave_arr = round(linspace(front_ave, rear_ave, total_point));

% SNR calculation noise added
snr_arr_noiseadd = nan(total_point-1, 1);
for i = 1:total_point-1
    Ascan_sub    = img_noiseadded(:,:,points_ave_arr(i):points_ave_arr(i+1));
    signal_power = mean(Ascan_sub.^2, "all");
    snr_arr_noiseadd(i)   = 10*log10(signal_power/noise_power);
end

close all;

% plot
figure,
plot(snr_arr);
hold on;
plot(snr_arr_noiseadd);

figure,
subplot(2,1,1);
image1 = img_3D(:,:,points_ave_arr(10));
imagesc(image1); colormap jet; colorbar;
subplot(2,1,2);
image2 = img_noiseadded(:,:,points_ave_arr(10));
imagesc(image2); colormap jet; colorbar;


figure,
subplot(2,1,1);
imagesc(img_3D(:,:,points_ave_arr(16))); colormap jet; colorbar;
subplot(2,1,2);
imagesc(img_noiseadded(:,:,points_ave_arr(16))); colormap jet; colorbar;

%% add noises to images directly
% SNR_arr = 15:25;
% dir_save = 'F:\Xiayang\results\Woven_samples\Dataset_learning';
% for snr = SNR_arr
%     for i = 1:total_point-1
%         img_batch = img_noiseadded(:,:,points_ave_arr(i)-1:points_ave_arr(i+1)+1);
%         img_batch = imnoise(img_batch,'gaussian', m, var_gauss);
%         save([dir_save  '\noisy_images\noisy_layer_' num2str(i) '_snr_' num2str(snr) '_' '.mat'], 'img_batch');
%     end
% end

%% add different snr and take 3 slices for each image

% save images as mat for training and testing
% noisy and desampled 
SNR_arr = 25;

dir_save = 'F:\Xiayang\results\Woven_samples\Dataset_learning';
for snr = SNR_arr
    sigPower       = sum(abs(img_3D).^2, 3) / lz; % linear
    reqSNR         = 10^(snr/10);
    noisePower     = sigPower/reqSNR;
    noise          = repmat(sqrt(noisePower), [1 1 lz]) .* randn(lx, ly, lz);
    img_noiseadded = img_3D + noise;
    for i = 1:total_point-1
        img_batch = img_noiseadded(:,:,points_ave_arr(i)-1:points_ave_arr(i+1)+1);
        save([dir_save  '\noisy_images\noisy_layer_' num2str(i) '_snr_' num2str(snr) '_' '.mat'], 'img_batch');
    end
end

%  origin
for i = 1:total_point-1
    img_batch = img_3D(:,:,points_ave_arr(i)-1:points_ave_arr(i+1)+1); % take 3 slices as average
    save([dir_save  '\origin_images_label\origin_layer_' num2str(i) '_' '.mat'], 'img_batch');
end

%% prepare dataset for LLC

dataset = img_3D(:,:,290:1080);
orthosliceViewer(dataset);

snr = 20;

sigPower       = sum(abs(img_3D).^2, 3) / lz; % linear
reqSNR         = 10^(snr/10);
noisePower     = sigPower/reqSNR;
lz_dataset =  size(dataset, 3);
noise          = repmat(sqrt(noisePower), [1 1 lz_dataset]) .* randn(lx, ly, size(dataset, 3));
dataset_noise = dataset + noise;

save('Ultrasound_dataset.mat', 'dataset', 'dataset_noise');
