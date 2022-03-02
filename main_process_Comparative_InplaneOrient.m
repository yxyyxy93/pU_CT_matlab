% read the file
clc;
close all;
fclose all;
clear;

%% Define the class
% read the preprocessed .mat
[Filename1,Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename = strcat(Pathname1, Filename1);
process = class_process_RFdata(filename);

% read the settings
[Filename1,Pathname1]=uigetfile({'*.xlsx'}, 'select the file');
settinig_file = strcat(Pathname1, Filename1);
process = process.loadSettings(settinig_file);

% read the data from the struct, or reset the dataset
process = process.read_origin_data;
% process = process.read_data; % read (reset) the dataset

% define the index to select Ascan 
x = 10;
y = 10;

%% EDA
process.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500];
window = [21 120 21 120 1 1500]; 
% window = [11 90 11 90 1 1500];  % for debug
process = process.cut_edges(window);
process.show_img_shape;
process.show_surfaces;

process.show_hilbert_Ascan(x, y);
process.show_Cscan(500, 'img_hil');

xslice =  50 / process.fx * 1e3;
yslice = 50 / process.fy * 1e3;
zslice = 1000 / process.fs * 1e6;
process.show_inph_3d('img_hil', xslice, yslice, zslice);

% filter and demo Scaleogram
f0 = linspace(1e6, 15e6, 100);
sigma = 0.8;
process.show_logGabor_Scaleogram(f0, sigma, x, y);

% filter the inph by logGabor filter
% filtered
f0 = 6.5e6;
sigma = 0.8;
process = process.Filter_logGabor(f0, sigma, 'img_hil');
process.show_inph_3d('img_hil_filter', xslice, yslice, zslice);
% calculate the unwraped inph 3D 
% smoothing scale 
sigma1 =[1, 1, 1];
process = process.show_unwraped_inph_3d('img_hil', xslice, yslice, zslice, sigma1);

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
process.show_Ascan_inam_peaks(20, 20, 50, 0.2);
% surface searching
process = process.find_front_amp(50, 0.2, 'img_hil');
process.show_surfaces;
% 
process = process.smooth_rear_I;
process.show_surfaces;

%% in-plane angles extract
% define a parallel C_scan_inam
center       = [40, 40];
radii        = 20;
dis_to_front = 0.85; % us
process      = process.perfect_rear_surface;
process      = process.define_parallel_inamCscan(dis_to_front, 'img', center, radii);

% rotate the Cscan image 
centers = 35:30:35 + 30*4;
ls = 5:5:25;
ws = 25:-5:5;
degree = 30;
for i = 1:length(centers)
    center = [centers(i), centers(i)];
    l = ls(i);
    w = ws(i);
    process = process.circle_rotate(l, w, center, degree);
end

% process = process.define_followUnwrapP_inamCscan(dis_to_front, 'img_hil');
% extract the 1D angular map along the depth for one point
radius             = 6:3:30;
theta              = 1:1:180;
x                  = 50;
y                  = 50;
distances_to_front = 0:0.01:3.7; % us
process.extract_local_ID_AlongDepth_RT('img_hil', radius, theta, x, y, distances_to_front);

% extract the inplane angles map from the parallel Cscan by RT
radius = 6:3:30;
theta = 1:1:180;
process = process.compute_orientation_by_RT(radius, theta);
% process = process.compute_orientation_by_RT_correct(radius, theta);
% process.show_ID_RT(51, 46, radius, theta); % y, x
% process.show_orientation_by_ID_RT(radius, theta);
% check diff. radius
for r_idx = 1:length(radius)
    process.show_orientation_by_ID_RT_oneradius(r_idx, radius, theta);
end

% extract 3D inplane angles by parallel Cscans, RT
distances_to_front = 0:0.01:3.7; % us
tic;
process = process.extract_local_orientation_RT_3D_parallel('img_hil', radius, theta, distances_to_front);
toc;
%
tic;
process = process.extract_local_orientation_modifiedRT_3D_parallel('img_hil', radius, theta, distances_to_front);
toc;

%%
% demonstrate the 3D orientaion
% show slice
xslice = 60 / process.fx * 1e3;
yslice = 60 / process.fy * 1e3;
zslice = [];
process.show_inplane_direction_3D(xslice, yslice, zslice);
process.show_inplane_direction_3D_modifyRT(xslice, yslice, zslice);

%% save and load part 'ID'
% extract the 1D angular map along the depth for each point by Gabor filter
wavelength = 6:3:30;
orientation = 1:1:180;
distances_to_front = 0:0.01:3.7; % us
foldname = ['F:\Xiayang\results\'...
    , Filename1(1:end - 5), '_ID'];
process.extract_local_ID_AlongDepth_LG('img_hil', wavelength, orientation, distances_to_front, foldname);
% read the .txt ID files
S = dir(fullfile(foldname, '*.txt'));
for ii = 1:numel(S)
    fileID = [foldname, '\', S(ii).name];
    M = dlmread(fileID, '\n');
    lgrow = M(1);
    lgcol = M(2);
    lentotal = lgrow * lgcol + 2;
    ID = NaN(floor(length(M) / lentotal), lgrow, lgcol);
    Angular_distribute = NaN(floor(length(M) / lentotal), lgcol);
    k = 1;
    while  lentotal * k + 1 <= length(M)
        ID_k = reshape(M(4 + lentotal * (k - 1) : lentotal * k + 1), [lgrow, lgcol]);
        ID(k, :, :) = ID_k;
        % sum up to form the 1D angular distribution
%         Angular_distribute(k, :) = max(ID_k, [], 1);
        Angular_distribute(k, :) = sum(ID_k, 1);
        k = k + 1;
    end
    % Derivative the 1D angular distribution
%     [Gx, Gy] = imgradientxy(Angular_distribute, 'sobel');
%     Gx_norm1 = normalize(Gx, 1);
%     Gx_norm2 = normalize(Gx, 2);
%     [Gx, Gy, Gz] = imgradientxyz(ID, 'sobel');
%     [Gmag,Gazimuth,Gelevation] = imgradient3(Gx,Gy,Gz);
    % 2D B mode display
%     figure('Name','sum_Gmag_2');
%     surf(squeeze(sum(Gmag, 2)));
%     figure('Name','sum_Gmag_2');
%     surf(squeeze(max(Gmag, [], 2)));
%     %
%     figure('Name','sum_Gz_2');
%     surf(squeeze(sum(Gz, 2)));  
%     figure('Name','max_Gz_2');
%     surf(squeeze(max(Gz, [], 2)));
%     figure('Name','sum_ID_2');
%     surf(squeeze(sum(ID, 2)));
%     figure('Name','max_ID_2');
%     surf(squeeze(max(ID, [], 2)));
%     colormap jet;
%     disp(ii);
end

%%
% Gabor fitler
wavelength  = 6:3:30;
orientation = 1:1:180;
% process = process.show_logGabor_filter('img_hil_filter', wavelength, orientation);
K           = 0.1;
process     = process.compute_logGabor_filter_withoutFig('img_hil', wavelength, orientation, K);
process.show_ID_logGabor_filter(95, 95, wavelength, orientation); % y, x
process.show_orientation_by_ID_allwl(wavelength, orientation);
process.show_orientation_by_ID_allwl_kurtosis(wavelength, orientation);
% check diff. wavelength
% for wl_idx = 1:length(wavelength)
%     process.show_orientation_by_ID_onewl(wl_idx, wavelength, orientation);
% end

% extract 3D inplane angles by parallel Cscans
distances_to_front = 0:0.01:3.7; % us
tic;
process = process.extract_local_orientation_3D_parallel_allwl('img_hil', wavelength, orientation, distances_to_front);
toc;

% show slice
xslice = 60 / process.fx * 1e3;
yslice = 60 / process.fy * 1e3;
zslice = [];
process.show_inplane_direction_3D_LG(xslice, yslice, zslice);

% extract 3D inplane angles by Cscans following the unwrapped phase. 
% not yet correct
% process = process.extract_local_orientation_3D_byuwp_allwl('img_hil', wavelength, orientation, distances_to_front);

%% save object
FolderName = "F:\Xiayang\results\";   % the destination folder;
savename = strcat(FolderName, Filename1(1:end-5), '_3D_orientation.mat');
process.save_local_orientaion_3D(savename);

%% read the ID .txt
[Filename1,Pathname1] = uigetfile({'*.txt'}, 'select the file');
filename = strcat(Pathname1, Filename1);
fprintf(filename,'%6f \r\n', A);

%% save all figures
FolderName = "F:\Xiayang\results";   % the destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  disp(fullfile(FolderName, FigName));
  set(0, 'CurrentFigure', FigHandle);
  saveas(gcf, strcat(FolderName, '\', FigName), 'bmp');
  saveas(gcf, strcat(FolderName, '\', FigName), 'fig');
end

%% save all figures
FolderName = "C:\Users\xiayang\OneDrive - UGent\matlab_work\process_RFcleaneddata\results";   % the destination folder
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