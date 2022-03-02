% read the file
% used to estimate thickness as well
% This is used for the test the orignal phase-derived interply track and the log-gabor fitelr
% for 1st paper parametric study

clc;
close all;
fclose all;
clear;

%%Define the class
% read the preprocessed .mat
[Filename1,Pathname1] = uigetfile({'*.mat'}, 'select the mat. file in ...\Yang_1stPaper_Para\input_data\exp');
filename              = strcat(Pathname1, Filename1);
process               = class_process_RFdata(filename);

% read the settings
[Filename1,Pathname1] = uigetfile({'*.xlsx'}, 'select the xlsx. file in ...\Yang_1stPaper_Para\input_data\exp');
settinig_file         = strcat(Pathname1, Filename1);
process               = process.loadSettings(settinig_file);

% read the data from the struct, or reset the dataset
process               = process.read_origin_data; % read (reset) the dataset

% define the index to select Ascan 
x                     = 15;
y                     = 15;

%% EDA
process.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500]; 
% window = [1 120 151 300 1 1500]; 
window         = [1 80 60 139 1 1500];  
x_step         = 1;
y_step         = 1;
process        = process.cut_edges(window, x_step, y_step);
process.show_img_shape;
% filename_fig   = filename;
% process.show_surfaces(filename_fig(1:end-4));

% process.show_hilbert_Ascan(x, y);
% process.show_Ascan_inam_peaks(150, 180, 50, 0.5);
% % surface calculation
% process        = process.find_front_amp(50, 0.5, 'img');
% process.show_surfaces(filename_fig(1:end-4));
% 
% xslice         = 50 / process.fx * 1e3;
% yslice         = 50 / process.fy * 1e3;
% zslice         = 1000 / process.fs * 1e6;
% process.show_inph_3d('img_hil', xslice, yslice, zslice);
% process.show_unwraped_inph_3d('img_hil', xslice, yslice, zslice);

% % filter and demo Scaleogram
% f0             = linspace(1e6, 10e6, 100);
% sigma          = 0.6;
% process.show_logGabor_Scaleogram(f0, sigma, x, y);
% %
% process.show_inph_3d('img_hil_filter', xslice, yslice, zslice);
% process.show_unwraped_inph_3d('img_hil_filter', xslice, yslice, zslice);

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
process.show_Ascan_inam_peaks(x, y, 300, 0.2, 'img_hil');
process.show_Ascan_inam_peaks(21, 30, 300, 0.2, 'img_hil');

% surface calculation
max_len      = 1500;
process      = process.find_front_amp(300, 0.2, 'img_hil', max_len, 1);
filename_fig = filename;
process.show_surfaces(filename_fig(1:end-4));
% %
% process      = process.smooth_rear_I;
% process.show_surfaces(filename_fig(1:end-4));
% %
% process      = process.recover_surface;

%% demonstrate A-scan and ply track by inph 
% interp_factor = 1;
% threshold     = 0.05;
% process.demo_Ascan_plytrack_inph(x, y, interp_factor, threshold);
% 
% % logGabor fitler added
% f0            = 6.5e6;
% sigma         = 0.7;
% process.demo_logGabor_plytrack_inph(x, y, f0, sigma, threshold);

%% demonstrate A-scan and ply track (noise added)
%
interp_factor = 1;
% process.demo_Ascan_plytrack(x, y, interp_factor);
process.demo_Ascan_plytrack_addnoise(x, y, interp_factor, 25);
% % surface echo removed, not good yet
% process.demo_removefront_plytrack(x, y);
% logGabor fitler added
f0            = 6.5e6;
sigma         = 0.7;
% process.demo_logGabor_plytrack(x, y, f0, sigma);
process.demo_logGabor_plytrack_addnoise(x, y, f0, sigma, 25)% surface echo removing + logGabor filter 
%
% process.demo_removefront_filter_plytrack(x, y, f0, sigma);

%% 3d ply track
% % add noise to img.
process    = process.addnoise(25);
% origin without filter
process    = process.track_interply('img_hil');
xslice     = 60/process.fx*1e3;
yslice     = 60/process.fy*1e3;
zslice     = [];
win        = 1:80;
% 
% process.show_track_interply(xslice, yslice, zslice);
win_x      = 1:80;
win_y      = 1:80; 
TOF_walls  = -mean(process.front_I-process.rear_I, 'all', 'omitnan');
TOF_oneply = TOF_walls/24; % 24 plies
process.show_oneinterply(5, 'nofilter', 3000, win_x, win_y, TOF_oneply);
process.show_oneinterply(23, 'nofilter', 3000, win_x, win_y, TOF_oneply);
process.show_B_scan_interply('x', 30, win, 'img_hil');

% filtered
process    = process.recover_surface;
f0         = 6.5e6;
sigma      = 0.7;
process    = process.Filter_logGabor(f0, sigma, 'img_hil');
process    = process.track_interply('img_hil_filter');
%
process.show_oneinterply(5,  'logGabor', 3000, win_x, win_y, TOF_oneply);
process.show_oneinterply(23, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process.show_B_scan_interply('x', 30, win, 'img_hil_filter');
%
figname      = 'logGabor';
process.thickness_estimation_v1(5.55, 24, win_x, win_y, figname);

%% save all figures
FolderName = "F:\Xiayang\results\compare_origin_filter";   % the destination folder
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