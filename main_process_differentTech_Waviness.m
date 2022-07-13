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
process3              = class_process_woven_RFdata(filename);

% read the settings
[Filename1,Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
settinig_file         = strcat(Pathname1, Filename1);
process3              = process3.loadSettings(settinig_file);

% read the data from the struct, or reset the dataset
process3              = process3.read_origin_data; % read (reset) the dataset

disp(process3.fx);
disp(process3.fy);

%% window cutting
% window                = [50 149 40 139 1 1500];
window                = [2 201 2 201 100 1400];
x_step                = 1;
y_step                = 1;
% window                = [2 300 2 480 1 1500];
% x_step                = 1;
% y_step                = 1;
process3              = process3.cut_edges(window, x_step, y_step);

process3.show_img_shape;

% shift the A scan in the time domain
process3 = process3.shift_A_scan(20);

%% surface calculation
% normal time window
PropertyName = 'img_hil';
% PropertyName = 'img_hil_filter';
delay        = 600;
max_len      = 1300;
flag         = 0;
front_I_max  = 800;
process3     = process3.find_damage_timewin_asignsurface(PropertyName, max_len, flag, delay, front_I_max, 1);
process3.show_surfaces(filename(1:end-5));

% process3.rear_I = process3.front_I + mean(TOF_walls);

% tackle the problem of walls determination
TOF_walls    = -mean(process3.front_I-process3.rear_I, 'all', 'omitnan');

% % shift the A scan in the time domain
% process3        = process3.shift_A_scan(200);

%% show A scan
% define the index to select Ascan 
x   = 50;
y   = 50;

[Ascan, t_space, fss] = process3.demo_Ascan(x, y, strcat('A_scan_', Filename1(1:end-5)));
   
%% ********** one-plane EDA
% define the C_scan_inam as well
% PropertyName = 'img_hil_filter';
PropertyName = 'img_hil';
ratio        = 3.5/24;
process3     = process3.define_parallel_inamCscan...
    (ratio, PropertyName);
%
imagename = 'C_scan_inam';
% curvelet denoising
sigma     = 5e-4;
process3  = process3.curvelet_denoise_inamCscan(sigma, imagename);
% % Radontransform
% theta         = 1:1:180;
% radiis        = 10;
% angle_compens = 0;
% process3      = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
% process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);

% 2d log-Gabor fitler
wavelength  = 4:4:16;
orientation = 0:2:180;
SFB         = 1; % [0.5 2.5]
SAR         = 0.5; % [0.23 0.92]

% imagename = 'C_scan_inam_denoise';
% process3  = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% process3.show_orientation_by_ID_allwl(wavelength, orientation, K);
imagename     = 'C_scan_inam';
process3      = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K             = 0.5e0;
angle_compens = -14;
tic;
process3.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);
toc;

%% 3D information Diagram
% PropertyName  = 'img_hil';
% 
% wavelength  = 8:2:16;
% orientation = 0:1:180; % degree
% SFB         = 1;       % [0.5 2.5]
% SAR         = 0.5;     % [0.23 0.92]
% z_theta     = 0;
% K           = 2;   
% z_range     = [1 1300];
% %
% tic;
% process3 = process3.extract_local_orientation_3DID...
%     (PropertyName, wavelength, orientation, z_theta, z_range, K, SFB, SAR);
% toc;
% 
% % % 3D ID plane
% % y_ori = 0:1:20;
% % z_ori = 0:1:20;
% % 
% % tic;
% % process3 = process3.extract_local_orientation_3DID_plane...
% %     (PropertyName, wavelength, y_ori, z_ori, z_range, K);
% % toc;
% 
% % show slice
% xslice        = 100 / process3.fx * 1e3;
% yslice        = 100 / process3.fy * 1e3;
% zslice        = [];
% mfsize        = [1 1 1];
% angle_compens = 14;
% process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);
% 
% process3.statistic_angular_distribution(angle_compens);

%% 'distance to front and rear' in-plane orientation extraction
% 2D log-Gabor filter
wavelength      = 4:4:16;
orientation     = 0:2:180;
ratio           = 0/24:0.3/24:24/24;
K               = 0.5e0;
% sigma              = 5e-4;
sigma           = 0;
tic;
process3 = process3.extract_local_orientation_3D_parallel_allwl(...
    PropertyName, wavelength, orientation, ratio, K, sigma);
toc;

% show slice
xslice        = 100 / process3.fx * 1e3;
yslice        = 100 / process3.fy * 1e3;
zslice        = [] / process3.fs; % all single type
mfsize        = [1 1 1];
angle_compens = -14;
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

process3.statistic_angular_distribution(angle_compens);

%% save dataset ***********
FolderName = "F:\Xiayang\results\DifferentTech_fiberorientation\";   % the destination folder
props = properties(process3);
for iprop = 1:length(props)
    thisprop = props{iprop};
    if(isnumeric(process3.(thisprop)))
        process3.(thisprop) = single(process3.(thisprop));
    end
end

save(strcat(FolderName, Filename1(1:end-5), '_50M.mat'), '-v7.3');

load(strcat(FolderName, '27072020_12h07m33s_50M.mat'));

% 27072020_14h22m13s_5M.mat
% 27072020_12h07m33s_50M

%% check A-scans - 2 times fft
close all;

process3.Cepstrum_deconv(50, 50);

%% save all figures
% "F:\Xiayang\results\Woven_samples\01032021_12h31m02s\f_5_omega07"
% "F:\Xiayang\results\Woven_samples\01032021_12h31m02s\f_10_omega07"
% !! don not forget to change the last folder path
FolderName = "F:\Xiayang\results\DifferentTech_fiberorientation\15MHz";   % the destination folder
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
