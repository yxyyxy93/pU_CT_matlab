classdef class_process_optimage
    % Addme
    % process the optical image
    % extract the local orientation
    
    properties
        img_origin
        img_gray
        % Radon transform
        RT_ID
        rear_mask
        % curvelet 
        Curvelet_ID
        Coef
        % Gabor filter - multi-resolution
        gaborMag
        gaborPha
        frow
        fcol
    end
    
    methods
        function obj = class_process_optimage(filename)
            % filename: path and name of the image
            obj.img_origin = imread(filename);
            if size(obj.img_origin,3)==3
                obj.img_origin = rgb2gray(obj.img_origin);
            else
                obj.img_origin = obj.img_origin;
            end
        end
        
        function show_img_shape(obj)
            % display the shapes of the properties
            props = properties(obj);
            for iprop = 1:length(props)
                thisprop = props{iprop};
                X_prop = [thisprop, ': ', num2str(size(obj.(thisprop)))];
                disp(X_prop);
            end
        end
        
        function show_figs(obj)
            % display the figs 
            figure('Name', 'original_fig');
            set(gcf, 'Position', [0, 0, 450, 350], 'color', 'white');
            ax = subplot(1, 1, 1);
            C_scan_inam_RT_pad = obj.img_origin;
            pcolor(ax, C_scan_inam_RT_pad);
            shading flat;
            set(gca, 'ydir', 'reverse');
            colormap gray;
            h  = colorbar;
            % caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} x (pixel)', 'fontsize', 16);
            ylabel('\fontname {times new roman} y (pixel)', 'fontsize', 16);
            % zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % **********************
            figure('Name', 'original_gray_fig');
            set(gcf, 'Position', [0, 0, 450, 350], 'color', 'white');
            ax = subplot(1, 1, 1);
            C_scan_inam_RT_pad = obj.img_gray;
            pcolor(ax, C_scan_inam_RT_pad);
            shading flat;
            set(gca, 'ydir', 'reverse');
            colormap gray;
            h  = colorbar;
            % caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} x (pixel)', 'fontsize', 16);
            ylabel('\fontname {times new roman} y (pixel)', 'fontsize', 16);
            % zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function obj = cut_edges(obj, window)
            % cut the edges of the dataset, and display the processed
            % figure
            % window: [x1, x2, y1, y2]
            x1 = window(1);
            x2 = window(2);
            y1 = window(3);
            y2 = window(4);
            obj.img_gray = obj.img_gray(x1: min(end, x2), y1: min(end, y2));
        end
        
        function obj = add_noise(obj, m, var_gauss)
            % adds Gaussian white noise with mean m and variance var_gauss.
            % m: mean of Gaussian white noise
            % var_gauss: variance of Gaussian white noise
            obj.img_gray = imnoise(obj.img_origin, 'gaussian', m, var_gauss);
%             obj.img_gray = imnoise(obj.img_origin, 'salt & pepper', var_gauss);
        end
        
        %%******************** radon transform
        function obj = compute_orientation_by_RT(obj, radius, theta, property_name)
            % extract the inplane orientation by Radon transform
            % radiis: the radius array for RT
            % theta: angle space
            % property_name: name of the image property
            % ***
            [lxc, lyc]         = size(obj.(property_name));
            ID                 = NaN(lxc, lyc, length(radius), length(theta));
            % ID_gyy = NaN(lxc, lyc, length(radiis), length(theta));
            % pad the C_scan_inam matrix
            r_max              = max(radius);
            C_scan_inam_RT_pad = obj.(property_name);
            % C_scan_inam_RT_pad = padarray(obj.img_gray, [r_max r_max], 'symmetric');
            % C_scan_inam_RT_pad = padarray(C_scan_inam_RT,[r_max r_max],mean(C_scan_inam_RT, 'all', 'omitnan'));
            % timer
            tic; 
            for r_idx = 1: length(radius)
                r = radius(r_idx);
                for i = 1 + r: lxc - r
                    for j = 1 + r: lyc - r
                        center                = [j, i];
                        % creat a round mask parameters
                        [~, anguler_1D, ~, ~] = fx_Radonto1Dangular_correct(C_scan_inam_RT_pad, center, r, theta);
                        ID(i, j, r_idx, :)    = anguler_1D;
                        % ID_gyy(i, j, r_idx, :) = anguler_1D_gyy;            
                    end
                    if mod(i, 10)==0
                        disp([num2str(i), '\', num2str(lxc)]);
                    end
                end
            end
            timeElapsed = toc;
            disp(['Compute_MultiResolution_RT:', num2str(timeElapsed)]);
            obj.RT_ID = ID;
        end
        
        function show_orientation_by_ID_RT(obj, radius, theta)
            % demonstrate the orientation extracted by Information diagram
            % of RT
            % use the inam image calcualated before
            % radius, theta: the same ratio for the multi-resolution RT
            % ***
            [lx, ly] = size(obj.img_gray);
            r_max    = max(radius);
            % remove the edges
            image_orientation = NaN(lx, ly);
            image_radius      = NaN(lx, ly);
            ID_sum            = zeros(length(radius), length(theta));
            for i = r_max+1: lx - r_max
                for j = r_max+1: ly - r_max
                    ID     = squeeze(obj.RT_ID(i-r_max, j-r_max, :, :));
                    maxval = max(ID);
                    if maxval <= mean(ID) + std(ID) % there is or is not significant fibrous content.
                        image_orientation(i, j) = NaN;
                    else
                        [~, indAng]             = max(ID);
                        image_orientation(i, j) = theta(indAng);
                    end
                    if length(radius)==1
                        [~, indAng] = max(ID);
                        indrad      = 1;
                    else
                        [~, I]           = max(ID(:));
                        [indrad, indAng] = ind2sub(size(ID),I);
                    end
                    image_orientation(i, j) = theta(indAng);
                    image_radius(i, j)      = radius(indrad);
                    % sum up the ID
                    ID_sum = ID_sum + ID;
                end
            end
            if ~isempty(obj.rear_mask)
                image_orientation(~obj.rear_mask) = NaN;
            end
            %
            if length(radius)~=1
                cf = figure();
                ax = subplot(1, 1, 1);
                surf(theta, radius, ID_sum);
                h = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
                set(ax, 'fontsize', 16);
                set(ax, 'linewidth', 1.5);
                xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
                ylabel('\fontname {times new roman} Wavelength (pixel)', 'fontsize', 16);
            end
            %
            figure('Name', ['orientation_image_RT_']);
            ax = subplot(1, 1, 1);
            X = (1: ly); % pixels
            Y = (1: lx); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            shading flat;
            colormap hsv;
            set(gca, 'YDir', 'reverse');
            caxis([0 180]);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
            %
            figure('Name', ['radius_image_RT_']);
            ax = subplot(1, 1, 1);
            X = (1: ly); % pixels
            Y = (1: lx); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_radius);
            shading flat;
            colormap jet;
            set(gca, 'YDir', 'reverse');
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Radius (pixel)');
        end
        
        function image_orientation = show_orientation_by_ID_RT_oneradius(obj, r_idx, radius, theta, pp_name)
            % demonstrate the orientation extracted by Information diagram
            % of RT
            % use the inam image calcualated before
            % r_idx: the index of the radius selected
            % radius, theta: the same ratio for the multi-resolution RT
            % pp_name: the name of the property of the image
            % ***
            [lx, ly]          = size(obj.(pp_name));
            r_max             = max(radius);
            % remove the edges
            image_orientation = NaN(lx, ly);
            r                 = radius(r_idx);
            for i = r+1: lx-r
                for j = r+1: ly-r
                    ID                      = squeeze(obj.RT_ID(i, j, r_idx, :));
                    %                     angular_1D = sum(ID, 1);
                    %                     [~, indAng] = max(angular_1D);
                    maxval = max(ID);
                    if maxval <= mean(ID) + std(ID) % there is or is not significant fibrous content.
                        image_orientation(i, j) = NaN;
                    else
                        [~, indAng]             = max(ID);
                        image_orientation(i, j) = theta(indAng);
                    end
 
                end
            end
            disp(r);
            if ~isempty(obj.rear_mask)
                image_orientation(~obj.rear_mask) = NaN;
            end
            %
            figure('Name', ['orientation_image_RT_' num2str(radius(r_idx))]);
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap hsv;
            caxis([0 180]);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \theta (\circ)');
        end
        
        function show_orientation_RT_allradius(obj, radius, theta)
            % demonstrate the orientation extracted by Information diagram
            % of RT
            % use the inam image calcualated before
            % r_idx: the index of the radius selected
            % radius, theta: the same ratio for the multi-resolution RT
            % ***
            [lx, ly]          = size(obj.img_gray);
            r_max             = max(radius);
            % remove the edges
            image_orientation = NaN(lx, ly);
            for i = r_max+1: lx-r_max
                for j = r_max+1: ly-r_max
                    ID                      = squeeze(obj.RT_ID(i-r_max, j-r_max, r_idx, :));
                    %                     angular_1D = sum(ID, 1);
                    %                     [~, indAng] = max(angular_1D);
                    [~, indAng]             = max(ID);
                    image_orientation(i, j) = theta(indAng);
                end
            end
            disp(radius(r_idx));
            if ~isempty(obj.rear_mask)
                image_orientation(~obj.rear_mask) = NaN;
            end
            %
            figure('Name', ['orientation_image_RT_' num2str(radius(r_idx))]);
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap hsv;
            caxis([0 180]);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
        end
        
        function show_1Dangular_RT(obj, x, y, r_idx, radius, theta)
            % display the information diagram of one pixel
            % y, x: the coordiantes of the pixel
            % r_idx: the index of the radius selected
            % radius, theta: the same ratio for the multi-resolution RT
            r_max              = max(radius);
            C_scan_inam_RT_pad = obj.img_gray;
            % get the local information concerning RT
            r      = radius(r_idx);
            center = [y, x];
            % creat a round mask parameters
            [R, anguler_1D, max_angle_I, xp] = fx_Radonto1Dangular_correct(C_scan_inam_RT_pad, center, r, theta);
            % ******  the disk-shaped area **********
            x1 = x - r;
            x2 = x + r;
            y1 = y - r;
            y2 = y + r;
            inam_C_scan_mask = C_scan_inam_RT_pad(x1: x2, y1: y2);
            [xx, yy]         = meshgrid(-r:r, -r:r);
            mask             = (hypot(xx, yy)<=r);
            inam_C_scan_mask = double(inam_C_scan_mask) .* mask;
            inam_C_scan_mask(mask==0) = NaN;
            figure('Name', ['disk-shaped_area' '_x_' num2str(y) '_y_' num2str(x)]);
            set(gcf, 'Position', [0, 0, 400, 300], 'color', 'white');
            ax = subplot(1, 1, 1);
            pcolor(ax, inam_C_scan_mask);
            shading flat;
            set(gca, 'ydir', 'reverse');
            hold on;
            colormap gray;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} x (pixel) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} y (pixel)', 'fontsize', 16);
            % zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % ******** the Radon transform formed from such projections
            figure('Name', ['Radon_transform_area' '_x_' num2str(y) '_y_' num2str(x)]);
            set(gcf, 'Position', [0, 0, 400, 300], 'color', 'white');
            ax = subplot(1, 1, 1);
            pcolor(ax, theta, xp, R);
            shading flat;
            set(gca, 'ydir', 'reverse');
            hold on;
            colormap gray;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} \theta (\circ)', 'fontsize', 16);
            ylabel('\fontname {times new roman} r (pixel)', 'fontsize', 16);
            % zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % ******** the 1D angular distribution by summing the absolute
            figure('Name', ['1D_angular_distribution' '_x_' num2str(y) '_y_' num2str(x)]);
            set(gcf, 'Position', [0, 0, 400, 300], 'color', 'white');
            ax = subplot(1, 1, 1);
            plot(ax, theta, anguler_1D, 'k-', 'linewidth', 2);
            hold on;
            scatter(ax, theta(max_angle_I), anguler_1D(max_angle_I), 'ro', 'filled', 'linewidth', 2)
            hold on;
            text(theta(max_angle_I), anguler_1D(max_angle_I), ...
                ['Maximum (dominant angle = ' num2str(theta(max_angle_I)) '^\circ)'], ...
                'Color', [1 0 0], 'FontSize', 16, 'Fontname', 'times new roman');
            xlim([theta(1) theta(end)]);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} \theta (\circ) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Mag. (arb.)', 'fontsize', 16);
            % zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % ***the original fig
            figure('Name', ['original_fig' '_x_' num2str(y) '_y_' num2str(x)]);
            set(gcf, 'Position', [0, 0, 400, 300], 'color', 'white');
            ax = subplot(1, 1, 1);
            C_scan_inam_RT_pad = obj.img_gray;
            pcolor(ax, C_scan_inam_RT_pad);
            shading flat;
            set(gca, 'ydir', 'reverse');
            hold on;
            scatter(y, x, 20, 'o', 'MarkerEdgeColor', 'red', 'LineWidth', 2);
            hold on;
            text(y+5, x, ['pixel at (' num2str(y) ',' num2str(x) ')'], ...
                'Color', 'red', 'FontSize', 16, 'Fontname', 'times new roman');
            colormap gray;
            h  = colorbar;
            % caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} x (pixel)', 'fontsize', 16);
            ylabel('\fontname {times new roman} y (pixel)', 'fontsize', 16);
            % zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        %%******************** Curvelet
        function obj = compute_curvelet(obj)
            % Forward curvelet transform
            disp('Take curvelet transform: fdct_usfft');
            tic; C     = fdct_wrapping(double(obj.img_gray), 0, 1); toc;
            img_temp_C = fdct_wrapping_dispcoef(C);
            % upsample the Coefficients by duplication
            % save as Information diagram
            [img_m, img_n] = size(obj.img_gray);
            ori_len_max    = length(C{end-1});
            C_temp         = NaN(img_m, img_n, length(C), ori_len_max);
            for i = 1:length(C)
                for j = 1:length(C{i})
                    curvelet_C        = abs(C{i}{j});  % absolute value in case of complex
                    % if the C is empty in last scale
                    if isempty(curvelet_C)
                        break;
                    end
                    len_ori           = ori_len_max / length(C{i});          
                    % using bulid-in 'imresize' function   
                    curvelet_C_resize = imresize(curvelet_C, [img_m img_n]);
%                     figure('Name', [num2str(i), '_' num2str(j)]), imagesc(curvelet_C_resize);
                    for idx_len_ori = 1:len_ori
                        C_temp(:, :, i, (j-1)*len_ori+idx_len_ori) ...
                            = curvelet_C_resize / 2^(-3*i/2); % theoretical scale
                        % In other words, one would need to compute on the order of 2^2j coefficients 
                        % per scale and angle as opposed to only about 2^3j/2 in the USFFT-based implementation.
                    end
%                     [m, n]     = size(curvelet_C);
%                     space_m    = round(img_m/m);
%                     space_n    = round(img_n/n);
%                     for im = 0:m-1
%                         for in = 0:n-1
%                             C_temp(space_m*im+1:space_m*(im+1), ...
%                                 space_n*in+1:space_n*(in+1), i, (j-1)*len_ori+1:j*len_ori) ...
%                                 = curvelet_C(im+1, in+1) / 2^(-3*i/2); % theoretical 
%                         end
%                     end 
% %                     figure; colormap gray; imagesc(abs(C{i}{j})); axis('image');
% %                     colorbar;
                end
            end   
            obj.Curvelet_ID = C_temp;
            obj.Coef        = C;
        end
        
        function show_curvelet_coef_onescale(obj, r_idx)
            % display the features of the curvelet coef.
            % ridx: the index for the selected scale
            [img_m, img_n]    = size(obj.img_gray); % get the size of the origin image
            C                 = obj.Coef;
            C_onescale        = C{r_idx};
            for i = 1:length(C_onescale)
                % using bulid-in 'imresize' function
                C_resize = imresize(C_onescale{i}, [img_m img_n]);
                cf = figure('Name', [num2str(i) '/' num2str(length(C_onescale))]);
                imagesc(abs(C_resize)); % absolute of complex
            end
        end
        
        function show_ID_curvelet(obj, x, y)
            % demonstrate the information diagram of curvelet transform 
            % and point in the original images
            % use the 4-D ID dataset calcualated before
            % x, y: index of the point in the image unit : pixel
            % ***
            ID             = squeeze(obj.Curvelet_ID(x, y, 1:end, :));
            angular_1D     = sum(ID, 1);
            % shift the angular_1D by 90 degrees (half length) in x axis
            %             angular_1D_shift = circshift(angular_1D, round(length(angular_1D)/2));
            %             angular_1D_sub = angular_1D - angular_1D_shift;
            [~, indAng]    = max(angular_1D);
            disp(180 / length(angular_1D) * indAng);
            % fillna
            C_scan_inam_LG = obj.img_gray;
            C_scan_inam_LG(isnan(C_scan_inam_LG)) = mean(C_scan_inam_LG, 'all', 'omitnan');
            %
            figure('Name', ['original image' ]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_LG, 2)); % default X = (1: size(C_scan_inam_LG, 2))  here
            Y = (1: size(C_scan_inam_LG, 1));
            imagesc(ax, X, Y, C_scan_inam_LG);
            hold on;
            scatter(ax, y, x, 30, 'r', 'filled'); % reversed x, y here
            colormap gray;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp.(arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman}', 'fontsize', 16);
            ylabel('\fontname {times new roman}', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % demonstrate the ID and 1D angular
            cf = figure();
            ax = subplot(1, 1, 1);
            imagesc(ID);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} orientaion ', 'fontsize', 16);
            ylabel('\fontname {times new roman} scale ', 'fontsize', 16);
        end
        
        function show_orientation_by_curveletID_onescale(obj, ridx)
            % demonstrate the orientation extracted by curvelet ID
            % use the 4D images calcualated before
            % ridx: the index for the selected scale
            % ***
            % kernel_size_half =
            [lx, ly]          = size(obj.img_gray);
            image_orientation = NaN(lx, ly);
            image_scale       = NaN(lx, ly);
            Curvelet_ID_4d    = obj.Curvelet_ID;
            for i = 1: lx
                for j = 1: ly
                    ID                      = squeeze(Curvelet_ID_4d(i, j, ridx, :));
                    [~, indAng]             = max(ID);
                    image_orientation(i, j) = mod((45 + indAng / size(Curvelet_ID_4d, 4) * 180), 180); % start from 45 degree
                end
            end
            %
            figure('Name', ['orientation_image']);
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap hsv;
            %             caxis([0 180]);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
        end
        
        function show_orientation_by_curveletID(obj)
            % demonstrate the orientation extracted by curvelet ID 
            % use the 4D images calcualated before
            % orientation: the orientation for Gabor filter, to
            % ***
            % kernel_size_half =
            [lx, ly]          = size(obj.img_gray);
            image_orientation = NaN(lx, ly);
            image_scale       = NaN(lx, ly);
            Curvelet_ID_4d    = obj.Curvelet_ID;
            for i = 1: lx
                for j = 1: ly
                    ID                      = squeeze(Curvelet_ID_4d(i, j, 1:end, :));
                    [~, I]                  = max(ID(:));
                    [indwl, indAng]         = ind2sub(size(ID), I);
                    image_orientation(i, j) = -mod((45 + indAng/size(Curvelet_ID_4d,4)*360), -180); % start from 45 degree
                    image_scale(i, j)       = indwl;
                end
            end
            %
            figure('Name', ['orientation_image']);
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap hsv;
%             caxis([0 180]);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
            %
            figure('Name', ['wavelength_image']);
            ax = subplot(1, 1, 1);
            X = (1: size(image_scale, 2)); % pixels
            Y = (1: size(image_scale, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_scale);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Wavelength (pixel)');
        end
        
        %%****************** Gabor filter ***********************8
        function obj = compute_orientation_Gabor_filter(obj, wavelength, orientation, SFB, SAR, property_name)
            % demonstrate and save the logGabor filtered images
            % use the C_scan_inam defined before
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % SFB, SAR: SpatialFrequencyBandwidth and SpatialAspectRatio
            % property_name: name of the image property
            % ***
            gaborArray = gabor(wavelength, orientation, ...
                'SpatialFrequencyBandwidth', SFB, 'SpatialAspectRatio', SAR);
            % fillna
            C_scan_inam_LG                        = obj.(property_name);
            NaN_part                              = isnan(C_scan_inam_LG);
            C_scan_inam_LG(isnan(C_scan_inam_LG)) = mean(C_scan_inam_LG, 'all', 'omitnan');
            wl_max                                = max(wavelength); % the maximum wavelength: used for cutting edge.
            C_scan_inam_RT_pad                    = C_scan_inam_LG; % no padding
            % % pad the C_scan_inam matrix
            % % C_scan_inam_RT_pad = padarray(C_scan_inam_LG, [round(wl_max / 2) round(wl_max / 2)], 'replicate');
            % C_scan_inam_RT_pad                    = padarray(C_scan_inam_LG, [round(wl_max / 2) round(wl_max / 2)], mean(C_scan_inam_LG, 'all', 'omitnan'));
            % %
%             figure;
%             ax = subplot(1, 1, 1);
%             X  = (1: size(C_scan_inam_RT_pad, 1));
%             Y  = (1: size(C_scan_inam_RT_pad, 2));
%             imagesc(ax, X, Y, C_scan_inam_RT_pad);
%             hold on;
%             h  = colorbar;
%             %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
%             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
%             set(ax, 'fontsize', 16);
%             set(ax, 'linewidth', 1.5);
%             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
%             ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
%             %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
%             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
%             set(gca, 'linewidth', 2);
            % Apply filters to C-scan image. **replication padding**
            [obj.gaborMag, obj.gaborPha] = imgaborfilt(...
                C_scan_inam_RT_pad, gaborArray);
            % normalized the magnitude.
            for i = 1:length(gaborArray)
                BW                    = gaborArray(i).SpatialFrequencyBandwidth;
                sigmaX                = gaborArray(i).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                sigmaY                = sigmaX ./ gaborArray(i).SpatialAspectRatio;
                gaborMag_i            = squeeze(obj.gaborMag(:, :, i));
                gaborMag_i            = gaborMag_i  / (2 * sigmaX * sigmaY * pi);
                gaborMag_i(NaN_part)  = NaN;
                obj.gaborMag(:, :, i) = gaborMag_i;
            end
            obj.frow = length(wavelength);
            obj.fcol = length(orientation);
        end
        
        function image_orientation = show_orientation_by_ID_allwl(obj, wavelength, orientation, K, property_name)
            % demonstrate the orientation extracted by Information diagram
            % use the logGabor filtered images calcualated before
            % orientation: the orientation for Gabor filter, to
            % determine the angle space
            % wavelength: the wavelength for Gabor filter, to determing the size of kernel for removing the edges
            % K:  the smoothing applied to the Gabor magnitude responses.
            % property_name: name of the image property
            % ***
            % kernel_size_half =
            [lx, ly]          = size(obj.(property_name));
            % remove the edges
            wl_max            = round(max(wavelength)/2);
            image_orientation = NaN(lx, ly);
            image_wavelength  = NaN(lx, ly);
            ID_sum            = zeros(length(wavelength), length(orientation));
            % Each Gabor magnitude image contains some local variations, 
            % even within well segmented regions of constant texture.
            % These local variations will throw off the segmentation.
            % We can compensate for these variations using simple Gaussian low-pass filtering to smooth the Gabor magnitude information.
            % We choose a sigma that is matched to the Gabor filter that extracted each feature.
            % We introduce a smoothing term K that controls how much smoothing is applied to the Gabor magnitude responses.
            gabormag = obj.gaborMag;
            g        = gabor(wavelength,orientation);  
            % Gaussian filter for the GaborMag
            if K~=0
                for i = 1:length(g)
                    sigma                    = 0.5*g(i).Wavelength;
                    gabormag_i               = gabormag(:, :, i);
                    invalid_part             = isnan(gabormag_i);
                    gabormag_i(invalid_part) = mean(gabormag_i, 'all', 'omitnan');
                    gabormag_i               = imgaussfilt(gabormag_i, K*sigma);
                    gabormag_i(invalid_part) = NaN;
                    gabormag(:,:,i)          = gabormag_i;
                end
            end
            for i = 1+wl_max: lx-wl_max
                for j = 1+wl_max: ly-wl_max
                    ID = squeeze(gabormag(i, j, :));
                    if isnan(ID(1, 1)) % invalid part has been assigned NaN
                        image_orientation(i, j) = NaN;
                        image_wavelength(i, j)  = NaN;
                    end
                    ID                      = reshape(ID, [obj.frow, obj.fcol]);
                    %                     ID_shift = circshift(ID, round(size(ID, 2) / 2), 2);
                    %                     ID_sub = ID - ID_shift;
                    [C, I]                  = max(ID(:));
                    if C <= mean(ID,'all') + std(ID, 0, 'all')
                        continue;
                    end
                    [indwl, indAng]         = ind2sub(size(ID),I);
                    image_orientation(i, j) = orientation(indAng);
                    image_wavelength(i, j)  = wavelength(indwl);
                    % sum up ID
                    ID_sum = ID_sum + ID;
                end
            end
            %             if ~isempty(obj.rear_mask)
            %                 image_orientation(~obj.rear_mask) = NaN;
            %             end
            %
%             cf = figure();
%             ax = subplot(1, 1, 1);
%             surf(orientation, wavelength, ID_sum);
%             h = colorbar;
%             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
%             set(ax, 'fontsize', 16, 'Fontname', 'times new Roman');
%             set(ax, 'linewidth', 1.5);
%             xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
%             ylabel('\fontname {times new roman} Wavelength (pixel)', 'fontsize', 16);
            %
            figure('Name', ['orientation_image_' property_name]);
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            set(gca, 'YDir', 'reverse');
            shading flat;
            xlabel('\fontname {times new roman} x (pixel)', 'fontsize', 16);
            ylabel('\fontname {times new roman} y (pixel)', 'fontsize', 16);
            colormap hsv;
            caxis([0 180]);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \theta (\circ)');
            set(ax, 'fontsize', 16, 'Fontname', 'times new Roman');
            set(ax, 'linewidth', 1.5);
            %
            figure('Name', ['wavelength_image_' property_name]);
            ax = subplot(1, 1, 1);
            X = (1: size(image_wavelength, 2)); % pixels
            Y = (1: size(image_wavelength, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_wavelength);
            set(gca, 'YDir', 'reverse');
            shading flat;
            xlabel('\fontname {times new roman} x (pixel)', 'fontsize', 16);
            ylabel('\fontname {times new roman} y (pixel)', 'fontsize', 16);
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \lambda (pixel)');
            set(ax, 'fontsize', 16, 'Fontname', 'times new Roman');
            set(ax, 'linewidth', 1.5);
        end
        
        function image_orientation = show_orientation_by_ID_onewl(obj, wl_idx, wavelength, orientation, K)
            % demonstrate the orientation extracted by one wavelength scale
            % orientation: the orientation for Gabor filter, to determine the angle space
            % wavelength: the wavelength for Gabor filter, to determing the size of kernel for removing the edges
            % K:  the smoothing applied to the Gabor magnitude responses.
            % SFBidx, SARidx: the indexes of the spatial frequency bandwidth and the spatial aspect ratio
            [lx, ly]          = size(obj.img_gray);
            % remove the edges
            wl_max            = round(max(wavelength));
            image_orientation = NaN(lx, ly);
            image_wavelength  = NaN(lx, ly);
            ID_sum            = zeros(length(wavelength), length(orientation));
            % Each Gabor magnitude image contains some local variations,
            % even within well segmented regions of constant texture.
            % These local variations will throw off the segmentation.
            % We can compensate for these variations using simple Gaussian low-pass filtering to smooth the Gabor magnitude information.
            % We choose a sigma that is matched to the Gabor filter that extracted each feature.
            % We introduce a smoothing term K that controls how much smoothing is applied to the Gabor magnitude responses.
            gabormag     = obj.gaborMag;
            g            = gabor(wavelength,orientation);
            % Gaussian filter for the GaborMag
            for i = 1:length(g)
                sigma                    = 0.5*g(i).Wavelength;
                gabormag_i               = gabormag(:, :, i);
                invalid_part             = isnan(gabormag_i);
                 if K~=0
                    gabormag_i(invalid_part) = mean(gabormag_i, 'all', 'omitnan');
                    gabormag_i               = imgaussfilt(gabormag_i, K*sigma);
                    gabormag_i(invalid_part) = NaN;
                end
                gabormag(:,:,i)          = gabormag_i;
            end
            for i = 1: lx - 2*wl_max
                for j = 1: ly - 2*wl_max
                    ID                                    = squeeze(gabormag(i, j, :));
                    if isnan(ID(1, 1)) % invalid part has been assigned NaN
                        image_orientation(i+wl_max, j+wl_max) = NaN;
                        image_wavelength(i+wl_max, j+wl_max)  = NaN;
                    end
                    ID                                    = reshape(ID, [obj.frow, obj.fcol]);
                    %                     ID_shift = circshift(ID, round(size(ID, 2) / 2), 2);
                    %                     ID_sub = ID - ID_shift;
                    [~, indAng]                           = max(ID(wl_idx, :));
                    image_orientation(i+wl_max, j+wl_max) = orientation(indAng);
                    image_wavelength(i+wl_max, j+wl_max)  = wavelength(wl_idx);
                    % sum up ID
                    ID_sum = ID_sum + ID;
                end
            end
            %             if ~isempty(obj.rear_mask)
            %                 image_orientation(~obj.rear_mask) = NaN;
            %             end
            %
%             cf = figure();
%             ax = subplot(1, 1, 1);
%             surf(orientation, wavelength, ID_sum);
%             h = colorbar;
%             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
%             set(ax, 'fontsize', 16);
%             set(ax, 'linewidth', 1.5);
%             xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
%             ylabel('\fontname {times new roman} Wavelength (pixel)', 'fontsize', 16);
            %
            figure('Name', ['orientation_image_wl_' num2str(wavelength(wl_idx))]);
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap hsv;
            caxis([0 180]);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
        end
        
        function show_ID_Gabor_filter(obj, x, y, wavelength, orientation)
            % display the information diagram of one pixel
            % y, x: the coordiantes of the pixel
            % wavelength, orientation: the parameters of the Gabor filter bank
            gabormag = obj.gaborMag;
            wl_max   = max(wavelength);
            % Gaussian filter for the GaborMag
            ID              = squeeze(gabormag(x, y, :));
            ID              = reshape(ID, [obj.frow, obj.fcol]);
            [~, I]          = max(ID(:));
            [indwl, indAng] = ind2sub(size(ID),I);
            figure('Name', ['ID_x_' num2str(x) '_y_' num2str(y)]);
            set(gcf, 'Position', [0, 0, 400, 300], 'color', 'white');
            ax = subplot(1, 1, 1);
            imagesc(ax, orientation, wavelength, ID);
            colormap gray;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');  
            hold on;
            scatter(ax, orientation(indAng), wavelength(indwl), 'ro', 'filled', 'linewidth', 2);
            hold on;
            text(orientation(indAng), wavelength(indwl), ...
                ['\leftarrow' 'Maxmimun (dominant angle = ' num2str(orientation(indAng)) '^\circ)'], ...
                'Color', [1 0 0], 'FontSize', 16, 'Fontname', 'times new roman');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} \theta (\circ)', 'fontsize', 16);
            ylabel('\fontname {times new roman} \lambda (pixel)', 'fontsize', 16);
            % ***the original fig
            figure('Name', ['original_fig' '_x_' num2str(y) '_y_' num2str(x)]);
            set(gcf, 'Position', [0, 0, 400, 300], 'color', 'white');
            ax = subplot(1, 1, 1);
            C_scan_inam_RT_pad = obj.img_origin;
            pcolor(ax, C_scan_inam_RT_pad);
            shading flat;
            set(gca, 'ydir', 'reverse');
            hold on;
            scatter(y, x, 20, 'o', 'MarkerEdgeColor', 'red', 'LineWidth', 2);
            hold on;
            text(y+5, x, ['\leftarrow' 'pixel at (' num2str(y) ',' num2str(x) ')'], ...
                'Color', 'red', 'FontSize', 16, 'Fontname', 'times new roman');
            colormap gray;
            h  = colorbar;
            % caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} x ', 'fontsize', 16);
            ylabel('\fontname {times new roman} y', 'fontsize', 16);
            % zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
    end
    
end

