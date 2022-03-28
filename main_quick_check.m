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

% read the settings from excel
[Filename1, Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
settinig_file          = strcat(Pathname1, Filename1);

% % read the data from the struct, or reset the dataset
% process3               = process3.loadSettings(settinig_file);

process3               = process3.read_origin_data; % read (reset) the dataset

%% cut window
window   = [51 300 151 400 1 1500];  
x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;

%% show A scan
% define the index to select Ascan 
x = 88;
y = 146;
% process3.demo_AS_3D_inclinfq(x, y);

process3.show_hilbert_Ascan(x, y);

%% % low-pass filter
bandpassFreq = 10e6;
process3     = process3.Filter_lowpass(bandpassFreq, 'img_hil');

%% surface search one signal
global min_pks 
min_pks      = 0.5;
PropertyName = 'img_hil';
% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y
clc;
close all;

MinPD   = 25;
MinPH   = 0.015;  % these 2 parameters need to be changed for surface estimation.
% surface calculation
max_len = 1400;
alpha   = 3e-3;
A_ratio = 0.6;

x = 88;
y = 146;
process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, x, y);

%% surface search all signals
% A_ratio      = 0.4; % for 0.3 m drop
process3     = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
close all;
process3.show_surfaces(filename_fig(1:end-5));

%% check 2 times fft
temp = process3.img;
temp = temp(130, 40, :);
temp = squeeze(temp);

figure, subplot(3, 1, 1);
plot(temp);
subplot(3, 1, 2);
temp_fft = fft(temp);
plot(abs(temp_fft));
subplot(3, 1, 3);
temp_fft2 = fft(abs(temp_fft));
plot(abs(temp_fft2));

process3 = process3.find_surface_2fft(MinPD, MinPH, PropertyName, max_len);


%% internal damage features
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