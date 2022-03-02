% read the file
clc;
close all;
fclose all;
clear;

%%
% 5 MHz
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

%% window cutting
% window                = [50 149 40 139 1 1500];
window                = [2 201 281 480 100 1400];
x_step                = 2;
y_step                = 2;
% window                = [2 300 2 480 1 1500];
% x_step                = 1;
% y_step                = 1;
process3              = process3.cut_edges(window, x_step, y_step);

process3.show_img_shape;

% shift the A scan in the time domain
process3 = process3.shift_A_scan(20);

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks

MinPD        = 300;
MinPH        = 0.1;  % these 2 parameters need to be changed for surface estimation.
PropertyName = 'img_hil';
process3.show_Ascan_inam_peaks(62, 53, MinPD, MinPH, PropertyName); % x, y

% surface calculation
max_len      = 1300;
threshold_d  = 1;
process3     = process3.find_front_amp(MinPD, MinPH, PropertyName, max_len, threshold_d);
filename_fig = filename;
process3.show_surfaces(filename_fig(1:end-5));
%
% process        = process.recover_surface;
process3     = process3.smooth_rear_I;
process3.show_surfaces(filename_fig(1:end-5));

% tackle the problem of walls determination
TOF_walls    = -mean(process3.front_I-process3.rear_I, 'all', 'omitnan');
% process3.rear_I = process3.front_I + mean(TOF_walls);
% process3.show_surfaces;

% % shift the A scan in the time domain
% process3        = process3.shift_A_scan(200);

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

% % log-Gabor filter
f0        = 6.5e6;
sigma     = 0.7;
process3  = process3.Filter_logGabor(f0, sigma, 'img_hil');

% % % low-pass filter
% bandpassFreq = 10e6;
% process3     = process3.Filter_lowpass(bandpassFreq, 'img_hil');

threshold    = 0.01;
% process   = process.track_interply('img_hil_filter');
process3  = process3.track_interply_inph(threshold, 'img_hil_filter', 25);

%%
% extract the scatter of the interply track with the markers indicating the
% layers.
% convert to 3d scattering points and initial the props for clustering
% process3 = process3.show_interply_track_3D;

%%
% % for debug
% threshold  = 0.05;
% x          = 80;
% y          = 140;
% process3.demo_logGabor_plytrack_inph_v2(x, y, f0, sigma, threshold);
% process3.show_logGabor_Scaleogram(3e6:0.1e6:8e6, 0.9, x, y);
        
B_type     = 'x';
index      = 80;
Bwin       = 1:100;
TOF_oneply = TOF_walls / 24; % 24 plies
process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil_filter');

% process.show_track_interply(xslice, yslice, zslice);
win_x     = 1:100;
win_y     = 1:100;
process3.show_oneinterply(05, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(15, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(22, 'logGabor', 3000, win_x, win_y, TOF_oneply);

std_mean_threshold = [1, 0.22]; % modify here 
process3.thickness_estimation(5.5, 24, win_x, win_y, std_mean_threshold, '5MHz_img_hil_filter');

% plot with uniformly spaced locaitons
nol = 25;
TOF = process3.track_show_B_scan_interply(...
    B_type, 80, Bwin, 'img_hil_filter', threshold, nol, std_mean_threshold);
% 
% % pos. err. calculation 
% % win_x      = 50:250;
% % win_y      = 50:250;
% layers_num = 25;
% process3   = process3.error_pos_calculate(win_x, win_y, layers_num);

% TOF    = process3.track_show_B_scan_interply(B_type, index, Bwin, 'img_hil_filter', threshold);
save TOF.mat TOF_walls;

%%
% 50 MHz
% Define the class
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

%% show A and B scan
fig_subname = '2-CEH105-24-p4_50MHz';
% define the index to select Ascan 
x           = 45;
y           = 45;
process.demo_Ascan(x, y, fig_subname);
% 
B_type      = 'y';
index       = 50;
Bwin        = 1:100;
process.demo_Bscan_interply(B_type, index, Bwin, fig_subname);

%%
window               = [2 501 2 501 100 1400];
x_step               = 5;
y_step               = 5;
process              = process.cut_edges(window, x_step, y_step);
process.show_img_shape;
% shift the A scan in the time domain
process        = process.shift_A_scan(200);

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
PropertyName = 'img_hil';
process.show_Ascan_inam_peaks(75, 46, 500, 0.02, PropertyName);
process.show_Ascan_inam_peaks(13, 30, 500, 0.02, PropertyName);
% surface calculation
max_len        = 1300;
process        = process.find_front_amp(500, 0.02, PropertyName, max_len);

load TOF.mat;
process.rear_I = process.front_I + mean(TOF_walls);
filename_fig   = filename;
process.show_surfaces(filename_fig(1:end-5));

%
% process = process.recover_surface;
% 
% process2.show_hilbert_Ascan(50, 55);

%% low-pass filter
bandpassFreq = 100e6;
process      = process.Filter_lowpass(bandpassFreq, 'img_hil');

%% ply track by inam 
%
property_name = 'img_hil_filter';
process       = process.track_interply_inam(property_name, 25);

B_type        = 'x';
index         = 50;
Bwin          = 1:100;
process.show_B_scan_interply(B_type, index, Bwin, property_name);

% plot with uniformly spaced locaitons
nol                = 25;
std_mean_threshold = [1, 0.22];
TOF                = process.track_show_B_scan_interply(B_type, index, Bwin, 'img', threshold, nol, std_mean_threshold);

% show ply track
win_x      = 1:100;
win_y      = 1:100;
TOF_oneply = TOF_walls / 24; % 24 plies
process.show_oneinterply(05, '50MHz', 3000, win_x, win_y, TOF_oneply);
process.show_oneinterply(15, '50MHz', 3000, win_x, win_y, TOF_oneply);
process.show_oneinterply(22, '50MHz', 3000, win_x, win_y, TOF_oneply);

% win_x      = 35:67;
% win_y      = 35:67;
std_mean_threshold = [1, 0.22];
process.thickness_estimation(5.5, 24, win_x, win_y, std_mean_threshold,'50MHz');

% % pos. err. calculation 
% layers_num = 25;
% process    = process.error_pos_calculate(win_x, win_y, layers_num, property_name);

%% **************************************
% 15 MHz
% Define the class
% read the preprocessed .mat
[Filename1,Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename              = strcat(Pathname1, Filename1);
process2              = class_process_RFdata(filename);

% read the settings
[Filename1,Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
settinig_file         = strcat(Pathname1, Filename1);
process2              = process2.loadSettings(settinig_file);

% read the data from the struct, or reset the dataset
process2              = process2.read_origin_data; % read (reset) the dataset
% 
% process2.show_hilbert_Ascan(101, 96);
% process2.show_hilbert_Ascan(50, 76);

%% window cutting
% process2.show_img_shape;
window   = [2 201 2 201 100 1400];
x_step   = 2;
y_step   = 2;
process2 = process2.cut_edges(window, x_step, y_step);
process2.show_img_shape;

%% read reference signal
[Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
filename               = strcat(Pathname1, Filename1);
x                      = 10;
y                      = 10;
process2               = process2.read_refer(filename, x, y);

% cut the ref. signal
win                    = [0.1e-6 2e-6]; % unit: s
process2               = process2.cut_reference_signal(win);

% calculate the SNR
win_noise              = [0.1e-6 0.5e-6];
win_signal             = [0.6e-6 1e-6]; % unit: s
snr                    = process2.calculate_SNR(win_signal, win_noise);

% align the signal
process2               = process2.align_refer_ascan(x, y); 

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
property_name = 'img';
process2.show_Ascan_inam_peaks(10, 22, 700, 0.03, property_name);
process2.show_Ascan_inam_peaks(78, 73, 700, 0.03, property_name);

% surface calculation
max_len      = 1300;
threshold_d  = 1;
MinPD        = 500;
MinPH        = 0.02;
process2     = process2.find_front_amp(MinPD, MinPH, PropertyName, max_len, threshold_d);
filename_fig = filename;
process2.show_surfaces(filename_fig(1:end-5));

process2.smooth_rear_I;
process2.show_surfaces(filename_fig(1:end-5));

% %
load TOF.mat;
% process2.rear_I = process2.front_I + mean(TOF_walls);
% % tackle the problem of walls determination
% process2.show_surfaces;

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
process2.demo_deconvolutions(x, y, q_factor, q_factor_AR, ...
    fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit, fig_subname);

% demo B scan
B_type      = 'x'; % scaning direction 'y'
index       = 50;
Bwin        = 1:100;
% process = process.deconv_ARextrapolation(Bwin, B_type, index, q_factor_AR, fft_padding, bw, bwr, k, ker_wid, fig_subname);

% process.demo_deconvolutions(x, y, q_factor, q_factor_AR, fft_padding, bw, k, lam, ker_wid, DownRate, Nit, fig_subname);
% 
process2     = process2.apply_deconvolutions_Bscan(Bwin, B_type, index, q_factor, q_factor_AR, fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit);

process2.demo_deconvolutions_Bscans(B_type, index, Bwin, fig_subname);

%%
% process2      = process2.apply_deconvolutions(q_factor, ...
%     q_factor_AR, fft_padding, bw, bwr, k, ker_wid);
process2      = process2.apply_deconvolutions_onlyWiener(q_factor, ker_wid);

% save('process2.mat','process2');

%% ply track by inam (deconvolve)
property_name = 'img_WienerDeconv';
process2      = process2.track_interply_inam(property_name, 25);

%
B_type        = 'x';
index         = 40;
Bwin          = 1:100;
TOF           = process2.track_show_B_scan_interply(B_type, index, Bwin, 'img_WienerDeconv', threshold, nol, std_mean_threshold);
% process2.show_B_scan_interply(B_type, index, Bwin, property_name);

% show ply track
win_x      = 1:100;
win_y      = 1:100;
TOF_oneply = TOF_walls / 24; % 24 plies
process2.show_oneinterply(05, property_name, 3000, win_x, win_y, TOF_oneply);
process2.show_oneinterply(15, property_name, 3000, win_x, win_y, TOF_oneply);
process2.show_oneinterply(22, property_name, 3000, win_x, win_y, TOF_oneply);

process2.thickness_estimation(5.5, 24, win_x, win_y, std_mean_threshold, '15MHz');

% % % % pos. err. calculation 
% nol        = 25;
% process2   = process2.error_pos_calculate(win_x, win_y, nol, property_name);

%% plot 3 by 3 images in one figure.
profile_layer               = cell(3, 3);
[~, ~, profile_layer{1, 1}] = process.show_oneinterply(05, property_name, 3000, win_x, win_y, TOF_oneply);
[~, ~, profile_layer{1, 2}] = process2.show_oneinterply(05, property_name, 3000, win_x, win_y, TOF_oneply);
[~, ~, profile_layer{1, 3}] = process3.show_oneinterply(05, property_name, 3000, win_x, win_y, TOF_oneply);
[~, ~, profile_layer{2, 1}] = process.show_oneinterply(15, property_name, 3000, win_x, win_y, TOF_oneply);
[~, ~, profile_layer{2, 2}] = process2.show_oneinterply(15, property_name, 3000, win_x, win_y, TOF_oneply);
[~, ~, profile_layer{2, 3}] = process3.show_oneinterply(15, property_name, 3000, win_x, win_y, TOF_oneply);
[~, ~, profile_layer{3, 1}] = process.show_oneinterply(22, property_name, 3000, win_x, win_y, TOF_oneply);
[~, ~, profile_layer{3, 2}] = process2.show_oneinterply(22, property_name, 3000, win_x, win_y, TOF_oneply);
[X, Y, profile_layer{3, 3}] = process3.show_oneinterply(22, property_name, 3000, win_x, win_y, TOF_oneply);
close all;

layers = {'5th interply', '15th interply', '22nd interply'};
techs  = {'50 MHz LP filtered signal', '15 MHz deconvolved signal', '5 MHz analytic-signal'}; 
cf     = fx_imagescs_tightplot(X, Y, layers, techs, profile_layer);

%% plot 3 images of thickness 'errorsbar' in one figure
win_x                      = 1:100;
win_y                      = 1:100;
errors                     = cell(1, 3);
stds                       = cell(1, 3);
std_mean_threshold         = [1 0.22];
[errors{1, 1}, stds{1, 1}] = process.thickness_estimation(5.5, 24, win_x, win_y, std_mean_threshold, '50MHz');
[errors{1, 2}, stds{1, 2}] = process2.thickness_estimation(5.5, 24, win_x, win_y, std_mean_threshold, '15MHz');
[errors{1, 3}, stds{1, 3}] = process3.thickness_estimation(5.5, 24, win_x, win_y, std_mean_threshold, '5MHz');
close all;

t_space                    = 1:1:24;
usthickness                = 5.52 / 24;
cf                         = fx_plot_estthick_tightplot...
    (t_space, techs, errors, stds, usthickness, std_mean_threshold);

%% plot 3 B-scan images of interplies tracking
xs           = cell(1, 3);
zs           = cell(1, 3);
Bscans_bg    = cell(1, 3);
Bscans_track = cell(1, 3);
front_lines  = cell(1, 3);
rear_lines   = cell(1, 3);
fss          = {process.fs, process2.fs, process3.fs};
fxs          = {process.fx, process2.fx, process3.fx};

B_type     = 'x';
index      = 80;
Bwin       = 1:100;
nol        = 25;

[~, xs{1}, zs{1}, Bscans_track{1}, Bscans_bg{1}, front_lines{1}, rear_lines{1}] = ...
    process.track_show_B_scan_interply(B_type, index, Bwin, 'img', threshold, nol, std_mean_threshold);
[~, xs{2}, zs{2}, Bscans_track{2}, Bscans_bg{2}, front_lines{2}, rear_lines{2}] = ...
    process2.track_show_B_scan_interply(B_type, index, Bwin, 'img_WienerDeconv', threshold, nol, std_mean_threshold);
[~, xs{3}, zs{3}, Bscans_track{3}, Bscans_bg{3}, front_lines{3}, rear_lines{3}] = ...
    process3.track_show_B_scan_interply(B_type, index, Bwin, 'img_hil_filter', threshold, nol, std_mean_threshold);
close all;

usthickness  = 5.5 / 24;
usTOF        = usthickness / 3000 * 1e3 * 2;
cf = fx_plot_bscans_tightplot(xs, zs, Bscans_bg, Bscans_track, ...
    front_lines, rear_lines, techs, fss, fxs, usTOF);

%% plot 3 A-scan recorded signal
x = 70; % best for 50 MHz
y = 70;

bandpassFreq   = 100e6;
process        = process.Filter_lowpass(bandpassFreq, 'img');
[ori_signals{1}, t_space, fss] = process.demo_Ascan(x, y, '2-CEH105-24-p4_50MHz');

x           = 40;
y           = 40;
% low-pass filter
bandpassFreq  = 30e6;
process2      = process2.Filter_lowpass(bandpassFreq, 'img');
[ori_signals{2}, ~] = process2.demo_Ascan(x, y, '2-CEH105-24-p4_15MHz');
bandpassFreq  = 10e6;
process3      = process3.Filter_lowpass(bandpassFreq, 'img');
[ori_signals{3}, ~] = process3.demo_Ascan(x, y, '2-CEH105-24-p4_5MHz');

close all;
win = [0.5 5];
cf  = fx_plot_ascans_tightplot(t_space, ori_signals, techs, fss, win);

%% demo errors of ply track 
% error bars
% PropertyName = 'img_hil';
% layers_num   = 25;
% threshold    = 0.15;
% process.demo_barErrors_phase(PropertyName, layers_num, threshold);
% 
% PropertyName = 'img_hil_filter';
% process.demo_barErrors_phase(PropertyName, layers_num, threshold);

% edge_value       = 0.225 / 0.9;
% edge_value_std   = 0.260 / 0.9;

edge_value       = 0.5;
edge_value_std   = 0.5;
fx_demo_best_results(process.errors_pos, process2.errors_pos, ...
    process3.errors_pos, edge_value, edge_value_std);

%% save all figures 1st paper
% FolderName = "F:\Xiayang\results\PlyTrack_inph";   % the destination folder
% FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   disp(fullfile(FolderName, FigName));
%   set(0, 'CurrentFigure', FigHandle); 
%   saveas(gcf, strcat(FolderName, '\', FigName), 'epsc');
%   saveas(gcf, strcat(FolderName, '\', FigName), 'pdf');
%   saveas(gcf, strcat(FolderName, '\', FigName), 'fig');
% end

%% save all figures 2nd paper
% FolderName  = "F:\Xiayang\results\DifferentTech_2";   % the destination folder
FolderName  = "F:\Xiayang\results\DifferentTech_2_revision";   % the destination folder
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

%% ************************** answers to questions ****************************




%% read reference signal
[Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
filename               = strcat(Pathname1, Filename1);
x                      = 10;
y                      = 10;
process               = process.read_refer(filename, x, y);

% cut the ref. signal
win                    = [0.1e-6 2e-6]; % unit: s
process               = process.cut_reference_signal(win);

% calculate the SNR
win_noise              = [0.1e-6 0.5e-6];
win_signal             = [0.6e-6 1e-6]; % unit: s
snr                    = process.calculate_SNR(win_signal, win_noise);

% align the signal
process               = process.align_refer_ascan(x, y); 

% process2.show_surfaces;

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
fig_subname = '2-CEH105-24-p4_50MHz';
process.demo_deconvolutions(x, y, q_factor, q_factor_AR, ...
    fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit, fig_subname);

%% 
x = 10;
y = 10;
process2.show_hilbert_Ascan(10, 10);
