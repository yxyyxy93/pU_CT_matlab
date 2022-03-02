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
path_load='F:\Xiayang\fromAGH\20170406_0757\';

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

%% ***********
FolderName = "F:\Xiayang\results\AGH_data\";   % the destination folder
C = strsplit(path_load,'\');
save(strcat(FolderName, C{end-1}, '.mat'), '-v7.3');

load(strcat(FolderName, '20170406_0757.mat'));
% '01032021_12h31m02s.mat'
% '01032021_14h17m04s.mat'
% '01032021_15h41m01s.mat'

%% window cut
% window = [101 180 81 160 1 1500]; 
% window = [1 120 151 300 1 1500]; 
window   = [101 300 101 300 1 650];  
x_step   = 1;
y_step   = 1;
process  = process.cut_edges(window, x_step, y_step);
process.show_img_shape;

%% show C scan
% define the index to select Ascan 
close all;
for z = 100:5:150
    % z            = ;
    PropertyName = 'img_hil';
%     process = process.show_Cscan(z, PropertyName);
    process.check2dfft(z, PropertyName); % 2d spectrum
end

%% 2d slice and 2d fft
% 
z            = 100;
PropertyName = 'img_hil';
% 
% 2d fft filter
% LP
close all;

process.check2dfft(z, PropertyName);

%% normalize along time domain;

process.normalize_timeAxis

%%
% 3D viewer
sliceViewer(abs(process.img_hil));

%% one-plane EDA
% define the C_scan_inam as well
% PropertyName = 'img_WienerDeconv';
% PropertyName = 'img_hil_filter';
PropertyName  = 'img_hil';
angle_compens = 0;

imagename     = 'C_scan_inam';

% Radontransform
theta       = 1:1:180;
radiis      = 20;
process    = process.compute_orientation_by_RT_correct(radiis, theta, imagename);

process.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);

[m_fiber_angle, std_fiber_angle, yhat, Idof_N] = process.calculate_m_std_fiber_angle_RT(radiis, theta, imagename, angle_compens);


%%
% 2d log-Gabor fitler
wavelength  = 8:4:40;
orientation = 1:1:180;
SFB         = 1; % [0.5 2.5]
SAR         = 0.5; % [0.23 0.92]
% imagename = 'C_scan_inam_denoise';
% process3  = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% process3.show_orientation_by_ID_allwl(wavelength, orientation, K);

tic;
PropertyName  = 'img_hil';
process    = process.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
toc;

% K: controls how much smoothing is applied to the Gabor magnitude responses.
tic;
K = 2e0;
process.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);
toc;

%% 
% monogenic signal

cw = 100;
tic;
process = process.compute_orientation_monogenicsignal(imagename, cw);
toc;


%% 'distance to front and rear' in-plane orientation extraction
% RT
% % PropertyName  = 'img_hil';

radius = 20;
theta  = 1:1:180; % degree
z      = 1:5:650;      
tic;
process       = process.extract_local_orientation_RT_3D_zaxis(...
    PropertyName, radius, theta, z);
toc;

% show slice
xslice        = 100 / process.fx * 1e3;
yslice        = 100 / process.fy * 1e3;
zslice        = []; % us
angle_compens = 0;
process.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);

% 2D log-Gabor filter
wavelength  = 8:4:40;
orientation = 1:1:180;
K           = 2e0;
% sigma              = 5e-4;
sigma       = 0;
tic;
process     = process.extract_local_orientation_3D_zaxis_allwl(...
    PropertyName, wavelength, orientation, z, K);
toc;

% show slice
xslice        = 100 / process.fx * 1e3;
yslice        = 100 / process.fy * 1e3;
zslice        = [];
mfsize        = [3 3 3];
angle_compens = 0;
process.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

%% ************** calculate the mean fiber angle and its standard deviation *************
% need the reference angle to calcualate the mean and std!

Stacking_sequence = [
    45 0 -45 -90 ...
    45 0 -45 -90 ...
    45 0 -45 -90 ...
    -90 -45 0 45 ...
    -90 -45 0 45 ...
    -90 -45 0 45];
m_fiber_angle_arr   = NaN(1, 24);
std_fiber_angle_arr = NaN(1, 24);
Idof_N_arr          = NaN(24, 180);
PropertyName_IFD    = 'Inplane_direction_3D_ID';
for p = 1:24 % 24 layers for default
    if p ==24
        ratio = (p-0.5) * 1/24;
    else
        ratio = (p-0.5) * 1/24;
    end
    ratio_next = (p-0.3) * 1/24;
    %     ref_angle  = Stacking_sequence(p);
    Idof_N_arr(p,:) = process3.calculate_angle_distrib(PropertyName, PropertyName_IFD, ratio, angle_compens);
    %     std_fiber_angle_arr(p) = std_fiber_angle_arr(p) / sum(Idof_Nr_arr(p,:));
    % normalize
    %     Idof_N_arr(p,:) = (Idof_N_arr(p,:)) / sum(Idof_N_arr(p,:));
    %     % log
    %     Idof_N_arr(p,:) = log(Idof_N_arr(p,:));
    disp([num2str(p) '/' '24']);
    %     if p==3 || p==11 || p==22
    if p==4 || p==5 || p==10 || p==15
        figure('Name', ['angle_distribution_' num2str(p) '_' PropertyName_IFD]);
        set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
        h1 = bar(-89-22:90-22, Idof_N_arr(p,:)/max(Idof_N_arr(p,:)));
        hold on;
        xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
        ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
        %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
        set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
        set(gca, 'linewidth', 2);
%         set(gca, 'YScale', 'log')
        %         legend([h1 h2], 'Angle distribution', 'Multiple Gaussian fitting', 'Interply track');
        %         ytickformat('percentage');
    end
end
%
maps             = cell(1,24);
for i = 1:24
    maps{i} = ['Ply' num2str(i) ':    ' num2str(Stacking_sequence(i)) '\circ'];
end
fx_creat_chartPlot(maps, Idof_N_arr, Stacking_sequence);

% m_fiber_angle_arr = m_fiber_angle_arr-90-23+angle_compens;
% 
m_fiber_angle_arr_RT   = NaN(1, 24);
std_fiber_angle_arr_RT = NaN(1, 24);
Idof_N_arr             = NaN(24, 180);
PropertyName_IFD       = 'Inplane_direction_3D_ID';
angle_compens          = -7;
for p = 1:24 % 24 layers for default
    if p ==24
        ratio = (p-0.5) * 1/24;
    else
        ratio = (p-0.5) * 1/24;
    end
    ratio_next = (p-0.3) * 1/24;
    ref_angle  = Stacking_sequence(p);
    [m_fiber_angle_arr_RT(p), std_fiber_angle_arr_RT(p), yhat, Idof_N_arr(p,:)] = process3.calculate_m_std_fiber_angle( ...
            PropertyName, PropertyName_IFD, ratio, ratio_next, angle_compens);
    %     std_fiber_angle_arr_RT(p) = std_fiber_angle_arr_RT(p) / sum(Idof_N_arr(p,:));
    Idof_N_arr(p,:) = process3.calculate_angle_distrib(PropertyName, PropertyName_IFD, ratio, angle_compens);
    %
    %     % normalize
%     Idof_N_arr(p,:) = (Idof_N_arr(p,:)) / sum(Idof_N_arr(p,:));
%     % log
%     Idof_N_arr(p,:) = log(Idof_N_arr(p,:));
    disp([num2str(p) '/' '24']);
    if p==3 || p==11 || p==22
        figure('Name', ['angle_distribution_' num2str(p) '_' PropertyName_IFD]);
        set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
        h1 = bar(-89-22:90-22, Idof_N_arr(p,:)/max(Idof_N_arr(p,:)));
        xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
        ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
        %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
        set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
        set(gca, 'linewidth', 2);
%         set(gca, 'YScale', 'log')
        %         legend([h1 h2], 'Angle distribution', 'Multiple Gaussian fitting', 'Interply track');
        %         ytickformat('percentage');
    end
end
fx_creat_chartPlot(maps, Idof_N_arr, Stacking_sequence);

% m_fiber_angle_arr_RT = m_fiber_angle_arr_RT-90-23 + angle_compens;
% 
% T = array2table([m_fiber_angle_arr; m_fiber_angle_arr_RT; ...
%     std_fiber_angle_arr; std_fiber_angle_arr_RT]);
% writetable(T,'test.xlsx','Sheet',1);

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

