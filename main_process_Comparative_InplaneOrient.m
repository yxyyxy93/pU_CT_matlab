% read the file
clc;
close all;
fclose all;
clear;

filevector = {'F:\Xiayang\2-ceh105-24-p5\v390\20220705_2-ceh105-24-p5_v390_25db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\v324\20220825_2-ceh105-24-p5_v324_22db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\v313\20220825_2-ceh105-24-p5_v313_23db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\v320\20220825_2-ceh105-24-p5_v320_17db_PE.mat', ...
    ... 'F:\Xiayang\2-ceh105-24-p5\v306\20220705_2-ceh105-24-p5_v306_21db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\geh5m\20220825_2-ceh105-24-p5_h5m_8db_PE.mat'};

fsvector      = [250 250 250 250 250 250] * 1e6;
windowsvector = {...
    [234 633 126 525 1 1400], ...
    [56 455 87 486 1 1500], ...
    [56 455 87 486 1 1500], ...
    [56 455 87 486 1 1400], ...
    ... [22 421 95 494 1 1400], ...
    [33 432 81 480 1 1400]};

%%

for i = 1:1
    %%
    %     i = 5;
    filename = filevector{i};
    
    %% Define the class
    % read the preprocessed .mat
    % [Filename1,Pathname1] = uigetfile({'*.mat'}, 'select the file');
    % filename = strcat(Pathname1, Filename1);
    process = class_process_RFdata(filename);
    
    % % read the settings
    % [Filename1,Pathname1]=uigetfile({'*.xlsx'}, 'select the file');
    % settinig_file = strcat(Pathname1, Filename1);
    % process = process.loadSettings(settinig_file);
    
    % read the data from the struct, or reset the dataset
    process = process.read_origin_data;
    % process = process.read_data; % read (reset) the dataset
    
    process.fx = 1 / 0.05e-3;
    process.fy = 1 / 0.05e-3;
    % process.fs = 250e6; % 50 MHz 25 MHz 15 MHz
    % process.fs = 125e6; % 7.5 MHz 5 MHz 2.5 MHz
    process.fs = fsvector(i);
    
    % define the index to select Ascan
    x = 10;
    y = 10;
    
    %% show A scan
    %     x   = 150;
    %     y   = 150;
    %
    %     [Ascan, t_space, fss] = process.demo_Ascan(x, y, strcat('A_scan_', 'x'));
    
    %% window
    % window   = [234 633 126 525 1 1400]; % for 20220705_2-ceh105-24-p5_v390_25db_PE.mat
    % window   = [33 432 116 515 1 1400]; % for 20220705_2-ceh105-24-p5_v324_20db_PE.mat
    % window   = [54 453 110 509 1 1400]; % for 20220705_2-ceh105-24-p5_v313_26db_PE
    % window   = [22 421 95 494 1 1400]; % for 20220705_2-ceh105-24-p5_v320_17db_1_PE
    % window   = [22 421 95 494 1 1400]; % for 20220705_2-ceh105-24-p5_v306_21db_PE.mat
    % window   = [23 422 100 499 1 1400]; % for 20220705_2-ceh105-24-p5_h5m_5db_PE.mat
    window   = windowsvector{i};
    x_step   = 2;
    y_step   = 2;
    process = process.cut_edges(window, x_step, y_step);
    process.show_img_shape;
    
    %% autocorrelation
    %     PropertyName = 'img_hil';
    %     process = process.caculate_autocorr(PropertyName);
    %
    %     orthosliceViewer(abs(process.img_autoc), 'Colormap', jet);
    
    %% normal time window
    PropertyName = 'img_hil';
    % PropertyName = 'img_hil_filter';
    delay       = 700;
    max_len     = 1400;
    flag        = 0;
    front_I_max = 500;
    process     = process.find_damage_timewin_asignsurface(PropertyName, max_len, flag, delay, front_I_max, 1);
    %     process.show_surfaces('surface');
    
    %% align the front surfaces
    %     delay_I      = 100;
    %     % PropertyName = 'img_hil';
    %     % PropertyName = 'img';
    %     process     = process.align_front_surface(delay_I, PropertyName);
    %
    %     process.rear_I = process.front_I + 3.8e-6 * process.fs; % 3.8 us;
    % %     back surface
    %     process3     = process3.align_rear_surface(delay_I, PropertyName);
    
    %% save space
%     process.rear_I    = process.front_I + 3.72e-6 * process.fs; % 3.8 us;
%     process.show_surfaces('surface');
%     % averaging testing; one-plane EDA
%     %     define the C_scan_inam as well
%     PropertyName = 'img_WienerDeconv';
%     PropertyName = 'img_hil_filter';
%     PropertyName  = 'img_hil';
%     %     PropertyName  = 'img_autoc';
%     angle_compens = -15;
%     % ***  ply 1
%     % ratio = 0.5 / 24;
%     % ratio = (0:0.05:1) / 24;
%     % *** ply 22
%     ratio = 21.5 / 24;
%     %
%     process = process.define_parallel_inamCscan...
%         (ratio, PropertyName);
%     imagename = 'C_scan_inam';
%     % % Radontransform
%     %     theta       = 1:1:180;
%     %     radiis      = 10;
%     %     process    = process.compute_orientation_by_RT_correct(radiis, theta, imagename);
%     %     process.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);
%     % 2d log-Gabor fitler
%     wavelength  = (4:2:16) * 2.5;
%     orientation = 1:1:180;
%     SFB         = 1; % [0.5 2.5]
%     SAR         = 0.5; % [0.23 0.92]
%     % imagename = 'C_scan_inam_denoise';
%     % process3  = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
%     % process3.show_orientation_by_ID_allwl(wavelength, orientation, K);
%     process    = process.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
%     % K: controls how much smoothing is applied to the Gabor magnitude responses.
%     K = 0.5e0;
%     process.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);
% %     
    %% 'distance to front and rear' in-plane orientation extraction
    % 2D log-Gabor filter
    %     PropertyName  = 'img_autoc';
    PropertyName  = 'img_hil';
    %     % RT
    %     radius        = 10;
    %     theta         = 1:1:180; % degree
    %     ratio         = 0/24:0.1/24:24/24;
    %     sigma_denoise = 0;
    %     yslice        = 50 / process.fy * 1e3;
    %     zslice        = []; % us
    %     tic;
    %     process       = process.extract_local_orientation_RT_3D_parallel(...
    %         PropertyName, radius, theta, ratio, sigma_denoise);
    %     toc;
    %     % show slice
    %     xslice        = 50 / process.fx * 1e3;
    %     angle_compens = -14;
    %     process.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);
    
    % ID
    wavelength  = (4:2:16) * 2.5;
    orientation = 1:1:180;
    ratio       = 0/24:0.1/24:24/24;
    K           = 0.5e0;
    % sigma              = 5e-4;
    sigma       = 0;
    tic;
    process     = process.extract_local_orientation_3D_parallel_allwl(...
        PropertyName, wavelength, orientation, ratio, K, sigma);
    toc;
    
    % show slice
    xslice        = 10 / process.fx * 1e3;
    yslice        = 10 / process.fy * 1e3;
    zslice        = [];
    
    mfsize        = [1 1 1];
    angle_compens = -15;
    process.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);
    
    %% ************** calculate the mean fiber angle and its standard deviation *************
    %     % need the reference angle to calcualate the mean and std!
    process.statistic_angular_distribution(angle_compens);
    
    %% save data
    Filename1  = strsplit(filename, '\');
    Filename1  = Filename1{end};
    FolderName = "F:\Xiayang\results\DifferentTech_fiberorientation\comparison\";   % the destination folder
    save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');
    
    %% save all figures
    FolderName = "F:\Xiayang\results\DifferentTech_fiberorientation\comparison\";   % the destination folder
    FolderName = strcat(FolderName, Filename1(1:end-4));
    mkdir(FolderName);
    FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
    for iFig = 1:length(FigList)
        FigHandle = FigList(iFig);
        FigName   = get(FigHandle, 'Name');
        disp(fullfile(FolderName, FigName));
        set(0, 'CurrentFigure', FigHandle);
        saveas(gcf, strcat(FolderName, '\', FigName), 'bmp');
        saveas(gcf, strcat(FolderName, '\', FigName), 'fig');
    end
    close all
end


% %%  ************************************ visualization part
% close all;
% clc;
%
% % show A scan
% x   = 100;
% y   = 100;
% [Ascan, t_space, fss] = process.demo_Ascan(x, y, strcat('A_scan_', 'x'));
% xlim([0.5 5.5]);
%
% % show C-scan slice
% PropertyName = 'img_hil';
% for ratio = [2.5 6.5 10.5] / 24
%     [~, C_scan_inam_temp, ~] = process.define_parallel_inamCscan(ratio, PropertyName);
%     x  = (1:size(C_scan_inam_temp, 1))/process.fx*1e3;
%     y  = (1:size(C_scan_inam_temp, 2))/process.fy*1e3;
%     figure, pcolor(x, y, C_scan_inam_temp);
%     shading flat;
%     colormap("jet");
%     colorbar;
%     set(gca, 'fontsize', 16);
%     set(gca, 'linewidth', 1.5);
% %     clim([0 0.12]);
%     if ratio ~= 10.5 / 24
%         axis off;
%     end
% end
%
% %%
% % show 2D orientation slice
% PropertyName = 'Inplane_direction_3D_ID';
% wl = max(wavelength);
% for ratio = [2.5 6.5 10.5] / 24
%     [~, C_scan_inam_temp, ~] = process.define_parallel_inamCscan(ratio, PropertyName);
%     x  = (1:size(C_scan_inam_temp, 1))/process.fx*1e3;
%     y  = (1:size(C_scan_inam_temp, 2))/process.fy*1e3;
%     C_scan_inam_temp = C_scan_inam_temp + angle_compens;
%     % remove the edges
%     C_scan_inam_temp(1:max(wl)/2+1, :) = NaN;
%     C_scan_inam_temp(:, 1:max(wl)/2+1) = NaN;
%     C_scan_inam_temp(end-max(wl)/2:end, :) = NaN;
%     C_scan_inam_temp(:, end-max(wl)/2:end) = NaN;
%     figure, pcolor(x, y, C_scan_inam_temp);
%     shading flat;
%     colormap("hsv");
%     colorbar;
%     clim([-90 90]);
%     set(gca, 'fontsize', 16);
%     set(gca, 'linewidth', 1.5);
%     if ratio ~= 10.5 / 24
%         axis off;
%     end
%     med  = median(C_scan_inam_temp(max(wl)/2+2:end-max(wl)/2-1,...
%         max(wl)/2+2:end-max(wl)/2-1),'all');
%     medd = mad(C_scan_inam_temp(max(wl)/2+2:end-max(wl)/2-1,...
%         max(wl)/2+2:end-max(wl)/2-1), 1, 'all');
%     disp([num2str(med) '+' num2str(medd)]);
% end
%
% %% 3D display
%
% % show slice
% xslice        = 100 / process.fx * 1e3;
% yslice        = 100 / process.fy * 1e3;
% zslice        = [];
% mfsize        = [1 1 1];
% angle_compens = -25;
% process.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);
%
% zlim([1 5.5]);
%
% %% ************** calculate the mean fiber angle and its standard deviation *************
% %     % need the reference angle to calcualate the mean and std!
% process.statistic_angular_distribution(angle_compens);
%
% ylim([1 5.5]);
% clim([0 0.5]);
