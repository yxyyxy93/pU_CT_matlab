% read the file

clc;
close all;
fclose all;
clear;

%% ***********
FolderName = "C:\Users\xiayang\OneDrive - UGent\matlab_work\results\Impact_damaged_sample\";   % the destination folder
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '14122020_12h08m51s.mat'));
% '20220207_BP1_3_025m_V313_25db_PEmat.mat'

%% Define the class
% read the preprocessed .mat
[Filename1, Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename               = strcat(Pathname1, Filename1);
process3               = class_process_woven_RFdata(filename);

% input the settings 
process3.fs = 250e6; % MHz
process3.fx = 1 / 0.25e-3; % unit: 1 / m
process3.fy = 1 / 0.25e-3; % unit: 1 / m

% % read the settings from excel
% [Filename1, Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
% settinig_file          = strcat(Pathname1, Filename1);
% % % read the data from the struct, or reset the dataset
% process3     = process3.loadSettings(settinig_file);

process3     = process3.read_origin_data; % read (reset) the dataset

sampling_rate = 1;

%% cut window
window   = [1 321 1 320 1 1500];  
x_step   = 3;
y_step   = 3;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;

%% original results display
xslice = 78 / process3.fx * 1e3;
yslice = 72 / process3.fy * 1e3;
zslice = 600 / process3.fs * 1e3 * 3000/2; % mm

process3.show_Raw_filtered_3D(xslice, yslice, zslice)

%% time varying gain
process3 = process3.read_origin_data; % read (reset) the dataset
ratio    = 1500;
process3 = process3.time_varying_gain(ratio);

%% show A scan
% define the index to select Ascan 
x   = 55;
y   = 75;
% process3.demo_AS_3D_inclinfq(x, y);
process3.show_hilbert_Ascan(x, y);

noises_dBW = [-30 -20 -10];
seed       = 1;
process3.show_Ascan_addnoise(x, y, noises_dBW, seed);

%% show A scan - sampling rate
x   = 32;
y   = 68;

samling_rates = [250/50 250/25 250/20];
process3.show_Ascan_resample(x, y, samling_rates)

%% read reference signal
[Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
filename               = strcat(Pathname1, Filename1);
x                      = 10;
y                      = 10;
process3               = process3.read_refer(filename, x, y);

% cut the ref. signal
win                    = [0.1e-6 3.5e-6]; % unit: s
process3               = process3.cut_reference_signal(win);

% align the signal
process3               = process3.align_refer_ascan(x, y); 

%% show A scan - different center frequencies
x   = 55;
y   = 75;
Fcs = [4, 3, 2] * 1e6;
process3.show_Ascan_freqManipulation(x, y, Fcs)

%% add noise
% process3 = process3.read_origin_data; % read (reset) the dataset
% process3 = process3.cut_edges(window, x_step, y_step);
% snr      = 15;
% process3 = process3.addnoise(snr);
close all;
dBW      = -10;
seed     = 1;
process3 = process3.addnoise_bydBW(dBW);

%% downsample dataset
% process3 = process3.read_origin_data; % read (reset) the dataset
% process3 = process3.cut_edges(window, x_step, y_step);
% snr      = 15;
% process3 = process3.addnoise(snr);
close all;
sampling_rate = 250/20;
process3      = process3.data_downsample(sampling_rate);

%% frequency manipulation dataset
clc;
close all;
Fc            = 5e6;
process3      = process3.data_freqManipulation(Fc);
sampling_rate = 1; % for reset the sampling_rate

%% surface search one signal
global min_pks 
min_pks      = 0.5;
% PropertyName = 'img_hil_filter';
% PropertyName = 'img_hil';
PropertyName = 'img_hil_noise';
% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y
clc;
close all;

MinPD   = 25;
MinPW   = 0.01;  % these 2 parameters need to be changed for surface estimation.
% surface calculation
max_len = 1500/sampling_rate;
alpha   = 3e-3;
A_ratio = 0.8;

x = 80;
y = 90;

process3.find_front_amp_alpha_Ascan(MinPD, MinPW, PropertyName, max_len, alpha, A_ratio, x, y);

%% surface search all signals
% A_ratio      = 0.4; % for 0.3 m drop
% process3     = process3.find_front_amp_alpha_maxfront(MinPD, MinPW, PropertyName, max_len, alpha, A_ratio);
process3     = process3.find_front_amp_alpha(MinPD, MinPW, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
close all;
 
process3.show_surfaces(filename_fig(1:end-5));

flag = 1; % using depth property % flag==1: using front and back surface to calculate
process3.damage_depth_imaging(flag);

%% deconvolution
clc
x = 140;
y = 140;

q_factor    = 1e-4;
q_factor_AR = 1e-4;
fft_padding = 2^11;
f1          = 3e6;
f2          = 7e6;
bw          = [f1 f2];
bwr         = -3:-0.5:-6;
k           = 15;
lam         = 0.01;
ker_wid     = 200;
Nit         = 45;
DownRate    = 1;
fig_subname = '2-CEH105-24-p4_5MHz';
process3.demo_deconvolutions(x, y, q_factor, q_factor_AR, ...
    fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit, fig_subname);

%%
% process3  = process3.apply_deconvolutions(q_factor, ...
%     q_factor_AR, fft_padding, bw, bwr, k, ker_wid);
process3      = process3.apply_deconvolutions_onlyWiener(q_factor, ker_wid);

% save('process2.mat','process2');

%% low pass filter
% PropertyName = 'img_hil';
% PropertyName = 'img_hil_noise';

% low-pass filter
Fc           = 5e6;
bandpassFreq = 2 * Fc;
process3     = process3.Filter_lowpass(bandpassFreq, PropertyName, sampling_rate);

%% normal time window
delay        = 60;
max_len      = 1500;
flag         = 1;
front_I_max  = 500;
% PropertyName = 'img_hil_filter';
% PropertyName = 'img_hil_noise';
PropertyName = 'img_hil';
% PropertyName = 'img_WienerDeconv';
process3.find_damage_timewin(PropertyName, max_len, flag, delay, front_I_max, sampling_rate);

%%
temp = process3.img;
temp = fft(temp, [], 3);
temp = abs(temp(:, :, 1:end/2));
temp = fft(temp, [], 3);
temp = abs(temp(:,:,1:end/2));
fx_Scrollable_3d_view(temp);

%% check 2 times fft
% close all;
clc;

x = 25 * process3.fx / 1e3;
y = 35 * process3.fy / 1e3;
% x = 55;
% y = 75;

% temp_img = real(process3.img_hil_noise);
temp_img = real(process3.img_hil);
% temp_img = real(process3.img_WienerDeconv);
% temp_img = real(process3.img_hil_filter);

temp   = temp_img(x, y, :);
temp   = squeeze(temp);
Fs_red = process3.fs/sampling_rate;
% 
L = length(temp);
L = 2*L;
% L = 2^nextpow2(L);

% DynamicGate = exp((1:L)/L);
% temp        = temp .* DynamicGate.';
% figure,
% % t_space = 1/Fs:1/Fs:L/Fs;
% plot(temp, 'LineWidth', 2);
% hold on;
% plot(abs(hilbert(temp)), 'LineWidth', 2);

temp_fft = fft(temp, L);
temp_fft = temp_fft(1:round(end/2));

temp_fft_abs = abs(temp_fft);

% temp_fft_abs = diff(temp_fft_abs, 1) ./ temp_fft_abs(2:end);
temp_fft_max = max(temp_fft_abs);

% *********** wavelet ************
wv = 'db3';
[c,l] = wavedec(temp_fft_abs, 4, wv);
% remove a4
c(1:l(1))    = 0 * c(1:l(1));
c(l(1):l(2)) = 0 * c(l(1):l(2));
c(l(2):l(3)) = 1 * c(l(2):l(3));
temp_rec = waverec(c, l, wv);
figure,
plot(temp_fft_abs, 'LineWidth', 2, 'DisplayName', 'spectrum');
hold on;
plot(temp_rec, 'LineWidth', 2, 'DisplayName', 'wavelet rec.');
hold on;
plot(diff(temp_fft_abs), 'LineWidth', 2, 'DisplayName', 'diff');
temp_fft_abs = temp_rec;
% ***********************

% temp_fft_abs(temp_fft_abs<temp_fft_max-3) = temp_fft_max-3;

% figure,
% plot(diff(temp_fft_abs), 'LineWidth', 2);

% % *********** band selection **********
% w = hann(round(10e6/Fs*L));
% temp_fft = temp_fft .* [w; zeros(length(temp_fft)-length(w), 1)];
% temp_fft(round(10e6/Fs*L):end) = 0;
% temp_fft(1:round(1e6/Fs*L))    = 0;
% temp_fft = temp_fft(round(1e6/Fs*L):round(8e6/Fs*L));
f_red        = Fs_red/L:Fs_red/L:Fs_red/2;
f_lo         = find(f_red<=1e6, 1, 'last' );f_up = find(f_red>=10e6, 1);
%             L_ori = 2^nextpow2(size(img_temp, 3)*sampling_rate);
L_ori        = L*sampling_rate;
temp_fft_abs = temp_fft_abs(f_lo:f_up);
% *************** find valley ****************
% findpeaks(-temp_fft_abs);
% figure, plot(valley_I);
% *************
temp_fft_abs = cat(1, temp_fft_abs, zeros(L_ori/2-(f_up-f_lo+1), 1));
figure,
plot(temp_fft_abs, 'LineWidth', 2);
% *********************

% f_space  = 0:Fs/L:Fs/2;
% temp_fft = cat(1, temp_fft, zeros(L*sampling_rate-length(temp_fft), 1));

temp_fft2_signal = fft(temp_fft_abs);
% t_space_2        = 0:2/Fs:L/Fs/2;
temp_fft2_signal = temp_fft2_signal(1:round(end/2));
figure,
plot(abs(temp_fft2_signal), 'LineWidth', 2);
% 
% % ************ reference signal
% A_scan_Ave    = squeeze(mean(temp_img(1:10,1:10,1:end/2), [1, 2]));
% % A_scan_Ave   = process3.refer_Ascan_aligned;
% A_scan_Ave   = real(A_scan_Ave);
% temp_fft_Ave = fft(A_scan_Ave, L);
% temp_fft_Ave = abs(temp_fft_Ave(1:round(end/2)));
% % *********** band selection **********
% L_ori    = L*sampling_rate;
% temp_fft_Ave = cat(1, temp_fft_Ave(f_lo:f_up), ones(L_ori/2-(f_up-f_lo+1), 1));
% 
% temp_fft2     = fft(abs(temp_fft_Ave), sampling_rate*L/2);
% temp_fft2_ave = abs(temp_fft2(1:round(end/2)));
% figure,
% plot(abs(temp_fft_Ave), 'LineWidth', 2);
% figure,
% plot((1:length(A_scan_Ave))/Fs, A_scan_Ave, 'LineWidth', 2);
% 
% % ******* divide
% temp_fft_abs = temp_fft_abs./temp_fft_Ave;
% figure,
% plot(temp_fft_abs, 'LineWidth', 2);
% figure,
% plot(diff(temp_fft_abs), 'LineWidth', 2);
% 
% temp_fft2_signal = fft(temp_fft_abs);
% % t_space_2        = 0:2/Fs:L/Fs/2;
% temp_fft2_signal = temp_fft2_signal(1:round(end/2));
% figure,
% plot(abs(temp_fft2_signal), 'LineWidth', 2);

% % substract
% figure,
% plot(abs(temp_fft2_signal)/max(abs(temp_fft2_signal)) ...
%     - abs(temp_fft2_ave)/max(abs(temp_fft2_ave)), 'LineWidth', 2);
% set(gca, 'fontsize', 16);
% set(gca, 'fontname', 'Times new roman');
% set(gca, 'linewidth', 1.5);
% 
%% 2 times fft determining the depth
PropertyName = 'img_hil';
max_len      = 1500;
flag_DG      = 0; % apply DynamicGate
delay        = 1; % delay for the search on the 2nd fft, unit: points
filter_flag  = 0;
passband     = [0.5e6 15e6];
process3     = process3.find_surface_2fft(PropertyName, max_len, flag_DG, delay, filter_flag, passband);

% filename_fig = filename;
% close all;
% process3.show_surfaces(filename_fig(1:end-5));

flag = 1; % using depth property % flag==1: using front and back surface to calculate
process3.damage_depth_imaging(flag);

%% 2 times fft determining the depth and imaging
% close all;
% PropertyName = 'img_hil';
PropertyName = 'img_hil_noise'; 
% PropertyName = 'img_hil_filter';
% PropertyName = 'img_WienerDeconv';

max_len     = 1500;
delay       = 40;
filter_flag = 0;
bwr         = -10; % dB

% process3.find_surface_2fft_imaging_diff(PropertyName, max_len, filter_flag, bwr, sampling_rate);
process3.find_surface_2fft_imaging_reference(PropertyName, max_len, filter_flag, bwr, sampling_rate);
% process3.find_surface_2fft_imaging_window(PropertyName, max_len, delay, filter_flag, bwr, sampling_rate);
% process3.find_surface_2fft_imaging_wavelet(PropertyName, filter_flag, sampling_rate);

%% check cepstrum
A = imread('C:\Users\xiayang\OneDrive - UGent\Pictures\Capture.PNG');
[IND, map] = rgb2ind(A, 2);
A_scan = nan(size(IND, 2), 1);
for i = 1:size(IND, 2)
    k = find(IND(:, i)==0);
    A_scan(i) = mean(k);
end
temp = A_scan(10:end-5);
temp = temp - mean(temp);

L = length(temp);

close all;

figure,

plot(temp, 'LineWidth', 2);
hold on;
plot(abs(hilbert(temp)), 'LineWidth', 2);

temp_fft = fft(temp);
figure, 
plot(abs(temp_fft), 'LineWidth', 2);

temp_fft = log(abs(temp_fft)) + 1i*angle(temp_fft);


temp_fft2_signal = ifft(temp_fft);

% % cut-off frequency
% cutoff_fp = 5;
% % temp_fft2_signal = temp_fft2_signal(cutoff_fp:end);
% temp_fft2_signal(1:cutoff_fp) = 0;
% temp_fft2_signal(end-cutoff_fp+1: end) = 0;

figure,
plot(real(temp_fft2_signal), 'LineWidth', 2);     

% inverse
temp_fft_invert = fft(temp_fft2_signal);
figure,
plot(abs(temp_fft_invert), 'LineWidth', 2);   
temp_fft_invert = exp(real(temp_fft_invert)) .* exp(1i*imag(temp_fft_invert));

temp_invert = ifft(temp_fft_invert);

figure,
plot(real(temp_invert), 'LineWidth', 2);     

%% check cepstrum 2 ************************
x = 55;
y = 55;
temp_img         = process3.img;
ascan            = squeeze(temp_img(x, y, :));
figure, plot(ascan);

ascan_fft        = fft(ascan);
ascan_fft        = ascan_fft(1:end/2);
figure, plot(abs(ascan_fft));
ascan_fft        = log(abs(ascan_fft)) + 1i*angle(ascan_fft);

ascan_fft2       = fft(ascan_fft);
figure, plot(real(ascan_fft2));

ascan_fft2(1:20) = 0;
% inverse
ascan_fft = fft(ascan_fft2(1:end/2));
ascan_fft = exp(real(ascan_fft)) .* exp(1i*imag(ascan_fft));
ascan     = ifft(ascan_fft, 'symmetric');
figure, plot(ascan);
        
%% time domain -reference determining the depth
PropertyName = 'img_hil';
max_len      = 1500;
flag_DG      = 0; % apply DynamicGate
delay        = 1; % delay for the search on the 2nd fft, unit: points
process3     = process3.find_surface_timedamin(PropertyName, max_len, flag_DG, delay);

filename_fig = filename;
close all;
process3.show_surfaces(filename_fig(1:end-5));

flag = 1; % using depth property % flag==1: using front and back surface to calculate
process3.damage_depth_imaging(flag);

%%
flag = 1; % using depth property
process3.damage_depth_imaging(flag);

%% freq-attenunation
close all;

x = 20;
y = 20;

temp_img = real(process3.img_hil);
temp     = temp_img(x, y, :);
temp     = squeeze(temp);
Fs       = process3.fs;
% 

front_fft = fft(temp(240: 390), 512);
rear_fft  = fft(temp(1160:1320), 512);

att = rear_fft ./ front_fft;
figure,
plot((1:512)*Fs/512, abs(att));


%% internal damage features
% % filtered
% f0        = 5e6;
% sigma     = 0.5;
% process3 = process3.Filter_logGabor(f0, sigma, 'img');

%
process3.damage_imaging;

%% save all figures  
FolderName  = "C:\Users\xiayang\OneDrive - UGent\matlab_work\results\quick_check";   % the destination folder
FigList     = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  % 
  FigName = [FigName '_flip']; % specific
  % 
  disp(fullfile(FolderName, FigName));
  set(0, 'CurrentFigure', FigHandle);
  saveas(gcf, strcat(FolderName, '\', FigName), 'epsc');
  saveas(gcf, strcat(FolderName, '\', FigName), 'pdf');
  saveas(gcf, strcat(FolderName, '\', FigName), 'fig');
end