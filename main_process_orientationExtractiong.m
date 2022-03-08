% read the file
clc;
close all;
fclose all;
clear;

%%
% 15 MHz Analytic-signal analysis
% Define the class
% read the preprocessed .mat
[Filename1, Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename               = strcat(Pathname1, Filename1);
process3               = class_process_RFdata(filename);

% read the settings
[Filename1, Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
settinig_file          = strcat(Pathname1, Filename1);

% read the data from the struct, or reset the dataset
process3               = process3.loadSettings(settinig_file);
process3               = process3.read_origin_data; % read (reset) the dataset


class_process_AGHdata

%% window cutting
% process2.show_img_shape;
window   = [51 350 1 300 1 1400];

window   = [21 220 21 220 1 1400];
x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;

% 
%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks
MinPD        = 800;
MinPH        = 0.02;  % these 2 parameters need to be changed for surface estimation.
% PropertyName = 'img_WienerDeconv'; bh
PropertyName = 'img_hil';
process3.show_Ascan_inam_peaks(153, 153, MinPD, MinPH, PropertyName); % x, y
% surface calculation
max_len                = 1400;
threshold_delamination = 0.3;
process3               = process3.find_front_amp(...
    MinPD, MinPH, PropertyName, max_len, threshold_delamination);
filename_fig           = filename;
% process3.show_surfaces(filename_fig(1:end-5));
%
process3.smooth_rear_I;
process3.show_surfaces(filename_fig(1:end-5));

% %
% load TOF.mat;
% process2.rear_I = process2.front_I + mean(TOF_walls);
% % tackle the problem of walls determination
% process2.show_surfaces;

%% original results display
% xslice = 10 / process3.fx * 1e3;
% yslice = 10 / process3.fy * 1e3;
% zslice = 2.4; % us
% process3.show_inaminph_3D(xslice, yslice, zslice);

%% averaging testing; one-plane EDA

% define the C_scan_inam as well
% PropertyName = 'img_WienerDeconv';
% PropertyName = 'img_hil_filter';
PropertyName  = 'img_hil';
angle_compens = -7;
% ***  ply 1
% ratio = 0.5 / 24;
% ratio = (0:0.05:1) / 24;
% *** ply 22
% ratio = 21.5 / 24;
ratio = (21.2:0.05:21.8) / 24;
%
process3    = process3.define_parallel_inamCscan...
    (ratio, PropertyName);
imagename     = 'C_scan_inam';
% % Radontransform
% theta       = 1:1:180;
% radiis      = 10;
% process3    = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
% process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);
% 2d log-Gabor fitler
wavelength  = 4:2:20;
orientation = 1:1:180;
SFB         = 1; % [0.5 2.5]
SAR         = 0.7; % [0.23 0.92]
% imagename = 'C_scan_inam_denoise';
% process3  = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% process3.show_orientation_by_ID_allwl(wavelength, orientation, K);
process3    = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K = 0.5e0;
process3.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);


%% 'distance to front and rear' in-plane orientation extraction
% RT
% % PropertyName  = 'img_hil';

radius        = 10;
theta         = 1:1:180; % degree
ratio         = 0/24:0.1/24:24/24;
sigma_denoise = 0;
tic;
process3      = process3.extract_local_orientation_RT_3D_parallel(...
    PropertyName, radius, theta, ratio, sigma_denoise);
toc;

% show slice
xslice        = 100 / process3.fx * 1e3;
yslice        = 100 / process3.fy * 1e3;
zslice        = []; % us
angle_compens = 0;
process3.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);

% 2D log-Gabor filter
wavelength         = 4:2:20;
orientation        = 1:1:180;
ratio              = 0/24:0.1/24:24/24;
K                  = 0.5e0;
% sigma              = 5e-4;
sigma              = 0;
tic;
process3           = process3.extract_local_orientation_3D_parallel_allwl(...
    PropertyName, wavelength, orientation, ratio, K, sigma);
toc;
% show slice
xslice        = 100 / process3.fx * 1e3;
yslice        = 100 / process3.fy * 1e3;
zslice        = [];
mfsize        = [1 1 1];
angle_compens = 0;
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

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

