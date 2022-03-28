% read the file
clc;
close all;
fclose all;
clear;

%%
% AGH University of Science and Technology
addpath('utils_AGH');

%%
% path to file
% path_load='F:\Xiayang\fromAGH\20170406_0757\';
path_load='C:\Users\xiayang\OneDrive - UGent\lab\experiment\fromAGH\20170406_0757\';

% read paramters of given scan
scanSet = readXMLscanner2([path_load 'info.xml']);

% number of files
n_files = length(dir([path_load '*.bin']));

%%
% filename               = strcat(Pathname1, Filename1);
process = class_process_AGHdata(path_load);

% input the settings 
process.fs = 200e6; % MHz
process.fx = 1 / 0.05e-3; % unit: 1 / m
process.fy = 1 / 0.05e-3; % unit: 1 / m

process = process.read_data_AGH(path_load, n_files);

%% window cut
% window = [101 180 81 160 1 1500]; 
% window = [1 120 151 300 1 1500]; 
window   = [1 1400 1 1400 1 650];  
x_step   = 4;
y_step   = 4;
process  = process.cut_edges(window, x_step, y_step);
process.show_img_shape;


%% normalize along time axis, remove DC component, and apply hilbert

process      = process.normalize_timeAxis_removeslope_hilbert;
% process      = process.removeslope_hilbert;

%% ***********
FolderName = "C:\Users\xiayang\OneDrive - UGent\matlab_work\results\AGH_data\";   % the destination folder
C = strsplit(path_load,'\');
save(strcat(FolderName, C{end-1}, '.mat'), '-v7.3');

load(strcat(FolderName, '20170406_0757.mat'));
% '01032021_12h31m02s.mat'
% '01032021_14h17m04s.mat'
% '01032021_15h41m01s.mat'

%% show A scan

% define the index to select Ascan 
x = 220;
y = 138;
process.show_hilbert_Ascan(x, y);

%% show B scan (filter)
B_type     = 'x';
index      = 140;
Bwin       = 1:350;
process.demo_Bscan_inst_filter(B_type, index, Bwin, 'img_hil');
    
%% show C scan
% define the index to select Ascan 
close all;
for z = 200:50:200
    % z            = ;
    PropertyName = 'img_hil';
%     process = process.show_Cscan(z, PropertyName);
    process.check2dfft(z, PropertyName); % 2d spectrum
end

%% surface searching
global min_pks % minimum amplitude for FWE
min_pks = 0.4;
            
x = 132;
y = 113;

PropertyName = 'img_hil';
% PropertyName = 'img_WienerDeconv';
MinPD        = 25;
MinPH        = 0.003;  % these 2 parameters need to be changed for surface estimation.
alpha        = 7e-3;
% A_ratio      = 0.9;
A_ratio      = 0.9;
max_len      = 650;
% max_len  = 1400;

% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y
% x = 111;
% y = 224;

% process3.show_Ascan_inam_peaks(x, y, MinPD, MinPH, PropertyName); % x, y
%
process.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, x, y)

%%
% surface calculation
% A_ratio      = 0.4; % for 0.3 m drop
process      = process.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio);
filename_fig = 'test';
process.show_surfaces(filename_fig(1:end-5));

% tackle the problem of walls determination
TOF_walls    = -mean(process.front_I-process.rear_I, 'all', 'omitnan');

%% amplitude rise method to determine the defect depth - spatial filter
clc;
close all;
PropertyName = 'img_hil';
zrange       = [450 500];
drop         = 3; % dB
sigma        = 1;

process = process.amplitude_rise_oneimage_spatialfilter(PropertyName, zrange, drop, sigma);

%% amplitude rise method to determine the defect depth 3D - spatial filter
PropertyName = 'img_hil';
zrange       = [1 650];
drop         = 3; % dB
sigma        = 1;
process = process.amplitude_rise_3D_spatialfilter(PropertyName, zrange, drop, sigma);
process.show_surfaces(filename_fig(1:end-5));

% process = process.rear_filter;

%% eda rear_I 
PropertyName = 'img_hil';
zrange       = [490 530];
process = process.search_rear_AS(PropertyName, zrange);

%% cluster algorithm
% K = 30;
% P = 2;
% nbins = 100;
% process = process.knn_search_rear_I(K, P, nbins);

%% correct rear surface by instantaneous phase
% 
% process = process.correct_rear_I_byinph;
% process.show_surfaces(filename_fig(1:end-5));

%% internal damage features
process.damage_imaging;

%% angular distribution analysis
clc;
PropertyName = 'img_hil';
zrange       = 1:1:650;
process = process.angulardistribution_2dfft(zrange, PropertyName);

%% 3d ply track
% 0.676 us - 2.13 us
% %  
% process = process.Filter_logGabor(f0, sigma, 'img');
% nol      = 49;
% % process3 = process3.track_interply_2ndharmonic('img_hil_filter', nol);
% process = process.track_interply('img_hil_filter');

% use 2nd-harmonic for first and the last interplies, 
% use fundamental resonance for another interplies

f0_1     = 12.8e6;
sigma0_1 = 0.8;
f0_2     = 6.3e6;
sigma0_2 = 0.7;
nol      = 26; % make it 25+1 in case that there is redundant accidently-tracked interply
process = process.track_interply_hybrid( ...
    'img', f0_1, sigma0_1, f0_2, sigma0_2, nol);


%%
B_type     = 'x';
index      = 140;
Bwin       = 1:350;
TOF_oneply = TOF_walls / 24; % 24 plies
process.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

%% structrul tensor 
% out-of-plane angles
f0    = 6.3e6;
sigma = 0.7;
process = process.Filter_logGabor(f0, sigma, 'img');

% smoothing scale
ds_rate      = 1;
sigma1       = [3, 3, 3];
% integration scale
sigma2       = [3, 3, 3];
PropertyName = 'img_hil_filter';
process      = process.structural_tensor(sigma1, sigma2, PropertyName, ds_rate);

% display the angles 3d slices
medf_kernel  = [3, 3, 3];
xslice       = 220 / process.fx * 1e3;
yslice       = 170 / process.fy * 1e3;
zslice       = [];
process.show_angles_ST_zangle(medf_kernel, xslice, yslice, zslice, ds_rate);

%
B_type     = 'x';
index      = 160;
Bwin       = 1:350;
process.show_angles_ST_Bscan(B_type, index, Bwin);

%
B_type     = 'y';
index      = 170;
Bwin       = 1:350;
process.show_angles_ST_Bscan(B_type, index, Bwin);

%% 
temp = process.c_p;
s = orthosliceViewer(temp, 'Colormap', jet, 'DisplayRange', [0.95 1]...
    , 'ScaleFactors', [1 1 0.5]);

temp_2D = sum(temp, 3);
figure, imagesc(temp_2D);

%% movie show


%% ***********
FolderName = "F:\Xiayang\results\image_processing\";   % the destination folder
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '13072020_16h54m52s.mat'));

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

