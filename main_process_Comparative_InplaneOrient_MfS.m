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
    
    delay_I     = 300;
    process     = process.align_front_surface(delay_I, PropertyName);
    process.show_surfaces('surface');
    
    %% Mumfordshah smoothing
    img_3D = abs(process.img_hil);
    img_3D = img_3D(:,:,300:round(mean(process.rear_I, 'all')));
    
    figure, imagesc(squeeze(img_3D(:, :, 500)));
    
    % set parameters
    alpha = 5;
    lambda = 0.1;
    eps = 1e-1;
    u_n = fx_mumfordshah(img_3D, alpha, lambda, eps);
    
    figure, imagesc(squeeze(u_n(:, :, 500)));
    
    %% ID
    % downsample
    wavelength  = (4:2:16) * 1;
    orientation = 1:1:180;
    K           = 0.5e0;
    
    SFB = 1;
    SAR = 0.5;
    gaborArray = gabor(wavelength, orientation, ...
        'SpatialFrequencyBandwidth', SFB, 'SpatialAspectRatio', SAR);
    
    u_n               = single(u_n(1:3:end, 1:3:end, :));
    Inplane_direction = NaN(size(u_n));
    
    for i = 1:2:size(u_n, 3)
        C_scan_inam_para_denoise = mean(u_n(:, :, i:i+1), 3);
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
        for ii = 1:length(gaborArray)
            sigma        = 0.5*gaborArray(ii).Wavelength;
            gabormag_i   = gaborMagnitude(:, :, ii);
            invalid_part = isnan(gabormag_i);
            if K~=0
                gabormag_i(invalid_part) = mean(gabormag_i, 'all', 'omitnan');
                gabormag_i               = imgaussfilt(gabormag_i, K*sigma);
                gabormag_i(invalid_part) = NaN;
            end
            gaborMagnitude(:,:,ii)          = gabormag_i;
        end
        % ********************** search for the max
        % *********** to accelerate, premise: average distributed orientation
        [~, Idx_1d] = max(gaborMagnitude, [], 3);
        indAng_2d   = floor(Idx_1d/length(wavelength));
        Ang         = orientation(1) + indAng_2d *(orientation(2) - orientation(1));
        Inplane_direction(:, :, i) = Ang;
        Inplane_direction(:, :, i+1) = Ang;
        disp(i);
    end
    
    %% ******************** 3D plot
    cf = figure('Name', ['3d_orientation_2Dfilter_', 'xslice', num2str(xslice(1))]);
    set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
    ax = subplot(1, 1, 1);
    y = (0: size(Inplane_direction, 1) - 1) / process.fx * 1e3;
    x = (0: size(Inplane_direction, 2) - 1) / process.fy * 1e3;
    z = (0: size(Inplane_direction, 3) - 1)/ process.fs * 1e3 * 3000/2; % depth axis
    [X, Y, Z] = meshgrid(x, y, z);
    %     % remove the edges
    %         inph_ex(1:max(wl)/2+1, :, :) = NaN;
    %         inph_ex(:, 1:max(wl)/2+1, :) = NaN;
    %         inph_ex(end-max(wl)/2:end, :, :) = NaN;
    %         inph_ex(:, end-max(wl)/2:end, :) = NaN;
    %
    Inplane_direction = mod(Inplane_direction + angle_compens, 180) - 90;
    h = slice(ax, X, Y, Z, Inplane_direction , xslice, yslice, zslice);
    hold on;
    set(h, 'EdgeColor', 'none');
    colormap hsv;
    h         = colorbar;
    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \theta (\circ)');
    set(h,'YTick', [-90 -45 0 45 90]); % set ticks
    caxis([-90 90]);
    set(ax, 'fontsize', 16);
    set(ax, 'linewidth', 1.5);
    xlabel('\fontname {times new roman} X (mm)', 'fontsize', 16);
    ylabel('\fontname {times new roman} Y (mm)', 'fontsize', 16);
    %                 zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
    zlabel('\fontname {times new roman} Time (\mum)', 'fontsize', 16); % time axis, unit: us
    set(gca, 'ZDir', 'reverse');
    
    xlim([x(1) x(end)]);
    ylim([y(1) y(end)]);
    %             legend([h3 h1 h2], 'Interply track', 'Front surface', 'Rear surface');
    %             ** Uncommnet for non-interply demp!
    xlim([0 x(end)]);
    ylim([0 y(end)]);
    %             zlim([0.8 7]);
    %             zticks([0 1 2 3 4 5 6]);
    view([15 65 40]);
    
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