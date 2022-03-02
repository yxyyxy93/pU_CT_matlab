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
process3               = class_process_IRTdata(filename);

% read the data from the struct, or reset the dataset
process3               = process3.read_origin_data_IRT; % read (reset) the dataset

%% original results display
% xslice = 20 / process3.fx * 1e3;
% yslice = 200 / process3.fy * 1e3;
% zslice = []; % us
% process3.show_inaminph_3D(xslice, yslice, zslice);
% 
% fx_Scrollable_3d_view(angle(process3.img_hil_filter));
% process3.show_hilbert_Ascan(100, 100);


%% log-gabor filter
% % filtered
% f0        = 7.5e6;
% sigma     = 0.7;
% process3 = process3.Filter_logGabor(f0, sigma, 'img');
% 
% % low-pass filter
% process3 = process3.Filter_lowpass(10e6, 'img_hil');
% 
% fx_Scrollable_3d_view(angle(process3.img_hil_filter));
% fx_Scrollable_3d_view(abs(process3.img_hil));

%% structrul tensor 

% smoothing scale
ds_rate      = 1;
sigma1       = [3, 3, 3];
% integration scale
sigma2       = [3, 3, 3];
PropertyName = 'img_hil';
process3     = process3.structural_tensor_IRT(sigma1, sigma2, PropertyName, ds_rate);

fx_Scrollable_3d_view(process3.img);

fx_Scrollable_3d_view(process3.c_p);
fx_Scrollable_3d_view(process3.c_l);

% % display the angles 3d slices
% medf_kernel  = [3, 3, 3];
% xslice       = 178 / process3.fx * 1e3;
% yslice       = 172 / process3.fx * 1e3;
% zslice       = [];
% process3.show_angles_ST_zangle(medf_kernel, xslice, yslice, zslice, ds_rate);
% %
% B_type     = 'x';
% index      = 200;
% Bwin       = 1:600;
% process3.show_angles_ST_Bscan(B_type, index, Bwin);
% %
% B_type     = 'y';
% index      = 300;
% Bwin       = 1:300;
% process3.show_angles_ST_Bscan(B_type, index, Bwin);

%% vedio of out-of-plane angles
B_type     = 'x';
Bwin       = 1:480;
process3.makeMovie_angles_ST_Bscan(B_type, Bwin);

%% ***********
FolderName = "F:\Xiayang\results\Woven_samples\";   % the destination folder
FolderName = "F:\Xiayang\results\Teijin_Wrinkling\";   % the destination folder
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '11032021_21h00m37s.mat'));
% '04042021_13h28m29s' 15 MHz
% '11032021_21h00m37s.mat' % 5 MHz

%% **************************************
% one-plane EDA
% 2D demo
% total_plies = 24;
% [angle, process] = process.show_one_Radontransform_vitualprofile(5, 0.5, 'img_hil_filter', center, radii, total_plies);

% % filtered
f0       = 15e6;
sigma    = 0.7; 
process3 = process3.Filter_logGabor(f0, sigma, 'img');

% define the C_scan_inam as well
% PropertyName = 'img_WienerDeconv';
PropertyName = 'img_hil_filter';
% PropertyName = 'img_hil';
ratio        = 8.5/24;
process3     = process3.define_parallel_inamCscan...
    (ratio, PropertyName);
%
imagename = 'C_scan_inam';
% curvelet denoising
sigma     = 5e-4;
process3  = process3.curvelet_denoise_inamCscan(sigma, imagename);

% Radontransform
theta         = 1:1:180;
radiis        = 10;
angle_compens = -7;
process3      = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);

% 2d log-Gabor fitler
wavelength  = 2:1:16;
orientation = 0:5:180;
SFB         = 1; % [0.5 2.5]
SAR         = 0.25; % [0.23 0.92]
% imagename = 'C_scan_inam_denoise';
% process3     = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% process3.show_orientation_by_ID_allwl(wavelength, orientation, K);
imagename     = 'C_scan_inam';
angle_compens = 0;
process3      = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K = 0.5e0;
process3.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);

%% 'distance to front and rear' in-plane orientation extraction
% RT
% PropertyName     = 'img_WienerDeconv'; % deconvolved signal
PropertyName  = 'img_hil_filter'; % orginal signal
% PropertyName  = 'img_hil';

radius        = 10;
theta         = 1:1:180; % degree
ratio         = 0/24:0.1/24:24/24;
sigma_denoise = 0;
tic;
process3      = process3.extract_local_orientation_RT_3D_parallel(...
    PropertyName, radius, theta, ratio, sigma_denoise);
toc;

% show slice
xslice        = 150 / process3.fx * 1e3;
yslice        = 150 / process3.fy * 1e3;
zslice        = []; % us
angle_compens = -7;
process3.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);

% 2D log-Gabor filter
wavelength         = 4:2:20;
orientation        = 0:5:180;
ratio              = 0/24:0.2/24:24/24;
K                  = 1e0;
% sigma              = 5e-4;
sigma              = 0;
tic;
process3           = process3.extract_local_orientation_3D_parallel_allwl(...
    PropertyName, wavelength, orientation, ratio, K, sigma);
toc;

% show slice
xslice        = 200 / process3.fx * 1e3;
yslice        = 200 / process3.fy * 1e3;
zslice        = [];
mfsize        = [1 1 1];
angle_compens = 0;
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

% 
B_type = 'x';
index  = 150;
Bwin   = 1:300;
process3.show_inplane_direction_3D_ID_Bscan(B_type, index, Bwin, [3 3 3], angle_compens);
%
B_type = 'y';
index  = 50;
Bwin   = 1:300;
process3.show_inplane_direction_3D_ID_Bscan(B_type, index, Bwin, [3 3 3], angle_compens);

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
    if p==3 || p==11 || p==22
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
PropertyName_IFD       = 'Inplane_direction_3D';
angle_compens          = -7;
for p = 1:24 % 24 layers for default
    if p ==24
        ratio = (p-0.5) * 1/24;
    else
        ratio = (p-0.5) * 1/24;
    end
    ratio_next = (p-0.3) * 1/24;
    ref_angle  = Stacking_sequence(p);
    %     [m_fiber_angle_arr_RT(p), std_fiber_angle_arr_RT(p), yhat, Idof_N_arr(p,:)] = process3.calculate_m_std_fiber_angle( ...
    %         PropertyName, PropertyName_IFD, ratio, ratio_next, angle_compens);
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

% %% save all figures 4th paper
% FolderName  = "F:\Xiayang\results\Impact_analysis";   % the destination folder
% FigList     = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   disp(fullfile(FolderName, FigName));
%   set(0, 'CurrentFigure', FigHandle);
%   saveas(gcf, strcat(FolderName, '\', FigName), 'epsc');
%   saveas(gcf, strcat(FolderName, '\', FigName), 'pdf');
%   saveas(gcf, strcat(FolderName, '\', FigName), 'fig');
% end

%% save to bin file

% matrix = process3.img;
% matrix = single(matrix);
% save('output.mat', 'matrix');
% 
% load('output.mat');
