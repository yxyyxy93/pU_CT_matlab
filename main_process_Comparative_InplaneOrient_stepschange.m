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

i = 3;
filename = filevector{i};

%% Define the class
% read the preprocessed .mat
% [Filename1,Pathname1] = uigetfile({'*.mat'}, 'select the file');
% filename = strcat(Pathname1, Filename1);
process = class_process_RFdata(filename);

%%
for step_ratio = [12 10 8 6 5 4 3]
    process = process.read_origin_data;
    % process = process.read_data; % read (reset) the dataset

    process.fx = 1 / (0.05e-3 * step_ratio);
    process.fy = 1 / (0.05e-3 * step_ratio);
    % process.fs = 250e6; % 50 MHz 25 MHz 15 MHz
    % process.fs = 125e6; % 7.5 MHz 5 MHz 2.5 MHz
    process.fs = fsvector(i);

    % define the index to select Ascan
    x = 10;
    y = 10;
    %% window
    % window   = [234 633 126 525 1 1400]; % for 20220705_2-ceh105-24-p5_v390_25db_PE.mat
    % window   = [33 432 116 515 1 1400]; % for 20220705_2-ceh105-24-p5_v324_20db_PE.mat
    % window   = [54 453 110 509 1 1400]; % for 20220705_2-ceh105-24-p5_v313_26db_PE
    % window   = [22 421 95 494 1 1400]; % for 20220705_2-ceh105-24-p5_v320_17db_1_PE
    % window   = [22 421 95 494 1 1400]; % for 20220705_2-ceh105-24-p5_v306_21db_PE.mat
    % window   = [23 422 100 499 1 1400]; % for 20220705_2-ceh105-24-p5_h5m_5db_PE.mat
    window   = windowsvector{i};
    x_step   = step_ratio;
    y_step   = step_ratio;
    process  = process.cut_edges(window, x_step, y_step);
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
    process.show_surfaces('surface');
    saveas(gcf, 'surface.bmp');
    close all;
    
    %% + Mumford Shah smoothing
        %     process.show_surfaces('surface');
    if i==1 % 50 MHz
        process.rear_I = process.front_I + 3.8e-6 * process.fs; % 3.8 us;
    end
    process.show_surfaces('surface');
    delay_I     = 300;
    process     = process.align_front_surface(delay_I, PropertyName);
    img_3D = abs(process.img_hil);
    img_3D = img_3D(:,:,300:round(mean(process.rear_I, 'all')));
    
    %% Mumfordshah smoothing
    if gpuDeviceCount
        img_3D = gpuArray(img_3D);
    end
    
    %     set parameters
    alpha = 2;
    lambda = 0.1;
    eps = 1e-1;
    u_n = fx_mumfordshah(img_3D, alpha, lambda, eps);
    img_3D = gather(img_3D);
    u_n = gather(u_n);
     
        %% ID
    % downsample
    wavelength  = (2:2:16) * 2;
    orientation = 1:1:180;
    K           = 0.5e0;
    
    SFB = 1;
    SAR = 0.5;
    gaborArray = gabor(wavelength, orientation, ...
        'SpatialFrequencyBandwidth', SFB, 'SpatialAspectRatio', SAR);
    
    u_n               = single(u_n(1:1:end, 1:1:end, :));
    Inplane_direction = NaN(size(u_n));
    
    for i = 1:1:size(u_n, 3)
        C_scan_inam_para_denoise = mean(u_n(:, :, i), 3);
        [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_para_denoise, gaborArray);
        gaborMagnitude      = single(gaborMagnitude);
        % normalized the magnitude.
        wl_Array       = repmat(wavelength', length(gaborArray)/length(wavelength), 1);
        BW             = SFB;
        sigmaX         = wl_Array / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
        sigmaY         = sigmaX ./ SAR;
        denominator    = permute(2 * sigmaX .* sigmaY * pi, [2, 3, 1]);
        gaborMagnitude = bsxfun(@rdivide, gaborMagnitude, denominator);
        % ********************************
        %         for ii = 1:length(gaborArray)
        %             sigma        = 0.5*gaborArray(ii).Wavelength;
        %             gabormag_i   = gaborMagnitude(:, :, ii);
        %             invalid_part = isnan(gabormag_i);
        %             if K~=0
        %                 gabormag_i(invalid_part) = mean(gabormag_i, 'all', 'omitnan');
        %                 gabormag_i               = imgaussfilt(gabormag_i, K*sigma);
        %                 gabormag_i(invalid_part) = NaN;
        %             end
        %             gaborMagnitude(:,:,ii)          = gabormag_i;
        %         end
        % ********************** search for the max
        % *********** to accelerate, premise: average distributed orientation
        [~, Idx_1d] = max(gaborMagnitude, [], 3);
        indAng_2d   = floor(Idx_1d/length(wavelength));
        Ang         = orientation(1) + indAng_2d *(orientation(2) - orientation(1));
        Inplane_direction(:, :, i) = Ang;
        disp(i);
    end
    
    %% release storate
    clear Ang;
    clear C_scan_inam_para_denoise;
    clear gabormag_i
    clear gaborMagnitude;
    clear gaborArray;
    clear Idx_1d;
    clear indAng_2d;
    clear invalid_part;
    clear img_3D;
    clear u_n;
    
    %% 'distance to front and rear' in-plane orientation extraction
%     % 2D log-Gabor filter
%     %     PropertyName  = 'img_autoc';
%     PropertyName  = 'img_hil';
%     %     % RT
%     %     radius        = 10;
%     %     theta         = 1:1:180; % degree
%     %     ratio         = 0/24:0.1/24:24/24;
%     %     sigma_denoise = 0;
%     %     yslice        = 50 / process.fy * 1e3;
%     %     zslice        = []; % us
%     %     tic;
%     %     process       = process.extract_local_orientation_RT_3D_parallel(...
%     %         PropertyName, radius, theta, ratio, sigma_denoise);
%     %     toc;
%     %     % show slice
%     %     xslice        = 50 / process.fx * 1e3;
%     %     angle_compens = -14;
%     %     process.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);
% 
%     % ID
%     wavelength  = (4:2:16) * 5 / step_ratio;
%     orientation = 1:1:180;
%     ratio       = 0/24:0.1/24:24/24;
%     K           = 1e0;
%     % sigma              = 5e-4;
%     sigma       = 0;
%     tic;
%     process     = process.extract_local_orientation_3D_parallel_allwl(...
%         PropertyName, wavelength, orientation, ratio, K, sigma);
%     toc;
% 
%     % show slice
%     xslice        = 200 / process.fx * 1e3 / step_ratio;
%     yslice        = 200 / process.fy * 1e3 / step_ratio;
%     zslice        = [];
%     mfsize        = [1 1 1];
%     angle_compens = -25;
%     process.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

%     %% ************** calculate the mean fiber angle and its standard deviation *************
%     %     % need the reference angle to calcualate the mean and std!
%     process.statistic_angular_distribution(angle_compens);

    %% save data
    Filename1  = strsplit(filename, '\');
    Filename1  = Filename1{end};
    FolderName = "F:\Xiayang\results\DifferentTech_fiberorientation\comparison\";   % the destination folder
    save(strcat(FolderName, Filename1(1:end-4), 'xyds_', num2str(step_ratio)), '-v7.3');

end

