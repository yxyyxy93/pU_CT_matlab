% read the file
clc;
close all;
fclose all;
clear;

filevector = {'F:\Xiayang\2-ceh105-24-p5\v390\20220705_2-ceh105-24-p5_v390_25db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\v324\20220825_2-ceh105-24-p5_v324_22db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\v313\20220825_2-ceh105-24-p5_v313_23db_PE.mat', ...
    ... 'F:\Xiayang\2-ceh105-24-p5\v320\20220825_2-ceh105-24-p5_v320_17db_PE.mat', ...
    ... 'F:\Xiayang\2-ceh105-24-p5\v306\20220705_2-ceh105-24-p5_v306_21db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\geh5m\20220825_2-ceh105-24-p5_h5m_8db_PE.mat'};

fsvector      = [250 250 250 250 250 250] * 1e6;
windowsvector = {...
    [234 633 126 525 1 1400], ...
    [56 455 87 486 1 1500], ...
    [56 455 87 486 1 1500], ...
    ... [56 455 87 486 1 1400], ...
    ... [22 421 95 494 1 1400], ...
    [33 432 81 480 1 1400]};

%%

for i = 1:4
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
    if i==1 % 50 MHz
        process.rear_I = process.front_I + 3.8e-6 * process.fs; % 3.8 us;
    end
    process.show_surfaces('surface');
    delay_I     = 300;
    process     = process.align_front_surface(delay_I, PropertyName);
    img_3D = abs(process.img_hil);
    img_3D = img_3D(:,:,300:round(mean(process.rear_I, 'all')));
    
    %     %% extract the 3D image by front-rear surface alignment
    %
    %     ratios = (0: 0.1: 24) / 24;
    %     [lx, ly, ~] = size(process.img_hil);
    %     img_3D = nan(lx, ly, length(ratios));
    %     for i = 1:numel(ratios)
    %         [~, C_scan_inam_para, ~] = process.define_parallel_inamCscan(ratios(i), PropertyName);
    %         img_3D(:, :, i) = C_scan_inam_para;
    %         disp(i);
    %     end
    %
    %     orthosliceViewer(img_3D, 'colormap', jet);
    
    %% Mumfordshah smoothing
    %     figure, imagesc(squeeze(img_3D(:, :, 500)));
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
    
    %         figure, imagesc(squeeze(u_n(:, :, 500)));
    
    %% test MumfoldShah
    close all;
    clc;
    
    lz = size(img_3D, 3);
    % 2D GF-ID
    wavelength  = (2:2:16) * 2;
    orientation = 1:1:180;
    SFB = 1;
    SAR = 0.5;
    gaborArray = gabor(wavelength, orientation, ...
        'SpatialFrequencyBandwidth', SFB, 'SpatialAspectRatio', SAR);
    angle_compens = -25;
    
    for ratio = [21.5] / 24
        %
        C_scan_inam_MS = u_n(:, :, round(ratio*lz));
        figure, imagesc(C_scan_inam_MS);
        [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_MS, gaborArray);
        gaborMagnitude      = single(gaborMagnitude);
        % normalized the magnitude.
        wl_Array       = repmat(wavelength', length(gaborArray)/length(wavelength), 1);
        BW             = SFB;
        sigmaX         = wl_Array / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
        sigmaY         = sigmaX ./ SAR;
        denominator    = permute(2 * sigmaX .* sigmaY * pi, [2, 3, 1]);
        gaborMagnitude = bsxfun(@rdivide, gaborMagnitude, denominator);
        [~, Idx_1d]    = max(gaborMagnitude, [], 3);
        indAng_2d      = floor(Idx_1d/length(wavelength));
        Ang_MS         = orientation(1) + indAng_2d *(orientation(2) - orientation(1));
        Ang_MS = mod(Ang_MS + angle_compens, 180) - 90;
        figure, imagesc(Ang_MS);
        colormap hsv;
        
        %
        C_scan_inam_temp = mean(img_3D(:, :, round(ratio * lz)-1:round(ratio * lz)+1), 3);
        [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_temp, gaborArray);
        gaborMagnitude      = single(gaborMagnitude);
        % normalized the magnitude.
        gaborMagnitude = bsxfun(@rdivide, gaborMagnitude, denominator);
        [~, Idx_1d] = max(gaborMagnitude, [], 3);
        indAng_2d   = floor(Idx_1d/length(wavelength));
        Ang         = orientation(1) + indAng_2d *(orientation(2) - orientation(1));
        Ang = mod(Ang + angle_compens, 180) - 90;
        figure, imagesc(C_scan_inam_temp);
        figure, imagesc(Ang);
        colormap hsv;
    end
    
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
    
    %% save data
    Filename1  = strsplit(filename, '\');
    Filename1  = Filename1{end};
    FolderName = "F:\Xiayang\results\DifferentTech_fiberorientation\comparison\";   % the destination folder
    save(strcat(FolderName, Filename1(1:end-4), 'MS', '.mat'), '-v7.3');
    
end
