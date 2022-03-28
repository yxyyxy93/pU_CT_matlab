classdef class_process_AGHdata < class_process_RFdata
    % sub class of class_process_RFdata
    properties
        defect_gt
        angular_distribution_fringerprint
    end
    
    methods
        function obj = class_process_AGHdata(filename)
            %
            obj = obj@class_process_RFdata(filename);
        end
        
        % ************************* display ***********************
        function demo_Bscan_inst_filter(obj, B_type, index, Bwin, PropertyName)
            % demonstrate the B_scan of the img
            % show inph_ex
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % Bwin: the indexes of the first and the last points in Bscan
            % PropertyName: the dataset to use
            % image
            img_temp = obj.(PropertyName);
            inph_ex  = angle(img_temp);
            % 3D Gaussian filter
            sigma1      = [3 3 3];
            inph_cos    = cos(inph_ex);
            %  inph_sin        = sin(inph_visual);
            inph_cos_filter = imgaussfilt3(inph_cos, sigma1, 'padding', 'replicate');
            %  inph_sin_filter = imgaussfilt3(inph_sin, sigma1, 'padding', 'replicate');
            inph_visual     = acos(inph_cos_filter);
            x = (1: Bwin(end)-Bwin(1)+1)  / obj.fx * 1e3;
            y = (1: Bwin(end)-Bwin(1)+1)  / obj.fy * 1e3;
            z = (1: size(inph_visual, 3)) / obj.fs * 1e6;
            %             ax = subplot(1, 1, 1);
            if (B_type == 'x')
                % B scan image
                img_temp = squeeze(img_temp(Bwin, index, :));
                img_temp = img_temp.';
                inph_ex  = squeeze(inph_ex(Bwin, index, :));
                inph_ex  = inph_ex.';
                inph_cos_filter  = squeeze(inph_cos_filter(Bwin, index, :));
                inph_cos_filter  = inph_cos_filter.';
                % surfaecs xslice
                if ~isnan(obj.front)
                    sur_f   = obj.front_I(Bwin, index)/ obj.fs * 1e6;
                    sur_r   = obj.rear_I(Bwin, index)/ obj.fs * 1e6;
                end
                y_index = x;
            elseif (B_type == 'y')
                % B scan image
                img_temp = squeeze(img_temp(index, Bwin, :));
                img_temp = img_temp.';
                inph_ex  = squeeze(inph_ex(index, Bwin, :));
                inph_ex  = inph_ex.';
                inph_cos_filter  = squeeze(inph_cos_filter(index, Bwin, :));
                inph_cos_filter  = inph_cos_filter.';
                % surfaecs yslice
                if ~isnan(obj.front)
                    sur_f   = obj.front_I(index, Bwin)/ obj.fs * 1e6;
                    sur_r   = obj.rear_I(index, Bwin)/ obj.fs * 1e6;
                end
                y_index = x;
                % choose the 'x' axis
                x       = y;
            end
            %             %remove the direct component
            %             B_scan = B_scan - mean(B_scan, 2);
            inam_ex        = abs(img_temp);
            inph_ex_unwarp = unwrap(inph_ex(end:-1:1, :), [], 1);  % reverse the array for correct differiatien
            infq_ex        = diff(inph_ex_unwarp, 1, 1) /2 / pi * obj.fs;
            infq_ex        = infq_ex(end:-1:1, :); % reverse bakc
            % ***** plot inam
            cf = figure('Name', ['Bscan' '_', B_type, '_', num2str(index), '_inam']);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            imagesc(x, z, inam_ex);
            hold on;
            % surfaecs slice
            if ~isnan(obj.front)
                h1 = scatter(y_index, sur_f, 2, 'red', 'filled', ...
                    'DisplayName','Front surface');
                hold on;
                h2 = scatter(y_index, sur_r, 2, 'magenta', 'filled', ...
                    'DisplayName','Rear surface');
                legend([h1 h2], 'Front surface', 'Rear surface');
                hold on;
            end
            % figure setup
            colormap(jet);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            cl = caxis;
            caxis(cl / 1.5);
            % ***** plot inph ************
            cf = figure('Name', ['Bscan' '_', B_type, '_', num2str(index), '_inph_cosfilter']);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            imagesc(x, z, inph_cos_filter);
            hold on;
            % surfaecs slice
            if ~isnan(obj.front)
                h1 = scatter(y_index, sur_f, 2, 'red', 'filled', ...
                    'DisplayName','Front surface');
                hold on;
                h2 = scatter(y_index, sur_r, 2, 'magenta', 'filled', ...
                    'DisplayName','Rear surface');
                legend([h1 h2], 'Front surface', 'Rear surface');
                hold on;
            end
            % figure setup
            colormap(gray);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Phase (rad.)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            % ***** plot infq
            cf = figure('Name', ['Bscan' '_', B_type, '_', num2str(index), '_infq']);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            imagesc(x, z(2:end), infq_ex/1e6);
            hold on;
            % surfaecs slice
            if ~isnan(obj.front)
                h1 = scatter(y_index, sur_f, 2, 'red', 'filled', ...
                    'DisplayName','Front surface');
                hold on;
                h2 = scatter(y_index, sur_r, 2, 'magenta', 'filled', ...
                    'DisplayName','Rear surface');
                legend([h1 h2], 'Front surface', 'Rear surface');
                hold on;
            end
            % figure setup
            colormap(jet);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Freq. (MHz)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            caxis([-20 30]);
        end
        
        % ************************* preprocess *************************
        
        function obj = read_data_AGH(obj, path_load, n_files)
            % loading raw data
            % read one file to get the data size
            [~, scanSet.noA, scanSet.noP] = loadBin([path_load '0avg.bin']);
            
            % allocatin  g space
            img_origin = single(zeros(n_files, scanSet.noA, scanSet.noP));
            
            % loop for loading
            for f_index = 0 : n_files-1
                % displaying progess
                clc
                disp(['Reading file ' num2str(f_index) ' of '  num2str(n_files) ])
                %load given file
                fn = [ num2str((f_index)) 'avg.bin'];
                [B, scanSet.noA, scanSet.noP] = loadBin([path_load fn]) ;
                %remove mean
                %     B = bsxfun(@minus,B,mean(B)); %disp('DC removal')
                %assign in 3D matrix
                img_origin(f_index+1, :, :) = B.';
            end
%             [lx, ly, lz] = size(img_origin);
%             img_hilbert = single(zeros(lx, ly, lz));
%             for i = 1:lx
%                 bscan     = squeeze(img_origin(i, :, :));
%                 bscan_hil = hilbert(bscan.'); %　hilbert for each colume at 2D matrix
%                 img_hilbert(i, :, :) = bscan_hil.'; %　hilbert for each colume at 2D matrix
%             
% %                 for j = 1:ly
% %                     img_hilbert(i, j, :) = hilbert(img_origin(i, j, :));
% %                 end
%                 clc
%                 fprintf('Hilbert progress: %0.2f%%\n',100*i/lx)
%             end
            obj.img = img_origin;
%             obj.img_hil = img_hilbert;
        end
        
        function obj = filter3D(obj, neighborhood)
            % filter the 3D　raw data
            img_origin = obj.img;
            img_origin = medfilter3(img_origin, neighborhood);
            obj.img = img_origin;
        end
        
        function obj = normalize_timeAxis_removeslope_hilbert(obj)
            img_temp    = obj.img;
            [lx, ly, ~] = size(img_temp);
            img_hilbert = nan(size(img_temp));
            for i = 1:lx
                for j = 1:ly
                    img_temp(i, j, :) = img_temp(i, j, :)/max(abs(img_temp(i, j, :)));
                    img_temp(i, j, :) = img_temp(i, j, :) - mean(img_temp(i, j, :));
                    img_hilbert(i, j, :) = hilbert(img_temp(i, j, :));
                end
                clc;
                fprintf('normalization progress: %0.2f%%\n',100*i/lx);
            end%;
            obj.img = img_temp;
            clear img_temp;
            obj.img_hil = img_hilbert;
            clear img_hilbert;
        end
        
        function obj = removeslope_hilbert(obj)
            %　no normalize along time axis
            img_temp    = obj.img;
            [lx, ly, ~] = size(img_temp);
            img_hilbert = nan(size(img_temp));
            for i = 1:lx
                for j = 1:ly
                    img_temp(i, j, :) = img_temp(i, j, :)/max(abs(img_temp(i, j, :)));
                    img_temp(i, j, :) = img_temp(i, j, :) - mean(img_temp(i, j, :));
                    img_hilbert(i, j, :) = hilbert(img_temp(i, j, :));
                end
                clc;
                fprintf('normalization progress: %0.2f%%\n',100*i/lx);
            end%;
            obj.img = img_temp;
            clear img_temp;
            obj.img_hil = img_hilbert;
            clear img_hilbert;
        end
        
        % *******************in-plane orientaion extraction plus *******
        
        function obj = compute_orientation_monogenicsignal(obj,imagename, cw)
            % extract the fiber orientation by monogenic signal
            % imagename: the 2D image name in the object
            % cw:  Centre-wavelengths in pixel units
            % ***
            I      = obj.(imagename);
            [Y, X] = size(I);
            % Now use these wavelengths to create a structure containing
            % frequency-domain filters to calculate the mnonogenic signal. We can
            % re-use this structure many times if we need for many images of the same
            % size and using the same wavelength. We can choose from a number of
            % different filter types, with log-Gabor ('lg') being the default. For lg
            % filters we can also choose the shape parameter (between 0 and 1), which
            % governs the bandwidth (0.41 gives a three-octave filter, 0.55 gives a two
            % octave filter)
            filtStruct = createMonogenicFilters(Y, X, cw, 'lg', 0.55);
            % Now we can use this structure to find the monogenic signal for the image
            [m1, m2, m3] = monogenicSignal(I, filtStruct);
            % The returned values are the three parts of the monogenic signal: m1 is
            % the even part, and m2 and m3 are the odd parts in the vertical and
            % horizontal directions respectively. Each array is Y x X x 1 x W, where
            % X and Y are the image dimensions and W is the number of wavelengths.
            % The filter responses to the filters of each scale are stacked along the
            % fourth dimension.
            %
            % (Alternatively one may pass a 3D volume to monogenicSignal, in which case
            % the 2D monogenic signal is found for each of the Z 2D planes independently and
            % returned as a set of Y x X x Z x W arrays)
            %
            % From here we can straightforwardly find many of the derived measures by
            % passing these three arrays
            %
            % Local energy (calculated on a per-scale basis)
            LE = localEnergy(m1, m2, m3);
            figure, imagesc(LE), colorbar;
            colormap(jet);
            % Local phase (calculated on a per-scale basis)
            LP = localPhase(m1, m2, m3);
            figure, imagesc(LP), colorbar;
            colormap(jet);
            %             Local orientation (calculated on a per-scale basis)
            %             Only need to pass the odd parts (m2,m3) as even part (m1) is irrelevant
            image_orientation = localOrientation(m2, m3);
            %             Display
            %             plot inplane orientation image
            cf = figure('Name', ['orientation_image_MS_' num2str(cw)]);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation/pi*90);
            shading flat;
            colormap hsv;
            set(gca, 'YDir', 'reverse');
            caxis([-90 90]);
            h = colorbar;
            h.Location = 'northoutside';
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \theta (\circ)');
            xlabel('\fontname {times new roman} x (pixel)', 'fontsize', 16);
            ylabel('\fontname {times new roman} y (pixel)', 'fontsize', 16);
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            % calcualate 1D distribution
            % minus the ref_angle to get the diff
            %             Inplane_direction_oneply = Inplane_direction_oneply - ref_angle;
            Idof_edge    = 0:180;
            Idof_N       = histcounts(image_orientation, Idof_edge);
            %             Idof_N       = circshift(Idof_N, 22); % shift 22 degree for fitting
            %************************* Gaussian fitting ********************
            centers      = [0 45 90 135];
            sigmas       = [8 8 8 8];
            amplitudes   = [1 1 1 1];
            %             centers      = [45 90 135];
            %             sigmas       = [8 8 8];
            %             amplitudes   = [1 1 1];
            numGaussians = 4;
            [amp, parameter, yhat] = fx_multi_Gaussianfitting(...
                numGaussians, Idof_edge(1:180), Idof_N, centers, sigmas, amplitudes);
            [C, amp_I]      = max(amp'.*parameter(2:2:end));
            m_fiber_angle   = parameter(amp_I*2-1);
            %             std_fiber_angle = amp(amp_I);
            %
            std_fiber_angle = std(image_orientation, 0, 'all','omitnan');
            disp(['mean:' num2str(m_fiber_angle) '. std:' num2str(std_fiber_angle)]);
            % ********** plot 1D distribution
            % normalize
            Idof_N = Idof_N / sum(Idof_N);
            % log
            figure('Name', ['angle_distribution_RT' '_' imagename]);
            set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
            h1 = bar(-89:90, Idof_N);
            xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
            %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
            set(gca, 'linewidth', 2);
            %         set(gca, 'YScale', 'log')
            %         legend([h1 h2], 'Angle distribution', 'Multiple Gaussian fitting', 'Interply track');
            %         ytickformat('percentage');
        end
        
        function [m_fiber_angle, std_fiber_angle, yhat, Idof_N] = calculate_m_std_fiber_angle_RT(obj, radius, theta, imagename, angle_compens)
            % calculate the mean fiber angle and its standard deviation
            % PropertyName: the property for defining the index
            % radius, theta: the same ratio for the multi-resolution RT
            [lx, ly]          = size(obj.(imagename));
            % remove the edges
            image_orientation = NaN(lx, ly);
            image_radius      = NaN(lx, ly);
            ID_sum            = zeros(length(radius), length(theta));
            r_max             = max(radius);
            for i = r_max+1: lx-r_max
                for j = r_max+1: ly-r_max
                    if isnan(obj.RT_ID(i, j))
                        continue;
                    else
                        ID = squeeze(obj.RT_ID(i, j, :, :));
                    end
                    if size(ID, 2)==1
                        [~, indAng] = max(ID);
                        indrad      = 1;
                    else
                        [~, I] = max(ID(:));
                        [indrad, indAng] = ind2sub(size(ID),I);
                    end
                    image_orientation(i, j) = theta(indAng) + angle_compens - 90;
                    image_radius(i, j)      = radius(indrad);
                    % sum up the ID
                    ID_sum                  = ID_sum + ID;
                end
            end
            %             Inplane_direction_oneply     = Inplane_direction;
            % minus the ref_angle to get the diff
            %             Inplane_direction_oneply = Inplane_direction_oneply - ref_angle;
            Idof_edge    = 0:180;
            Idof_N       = histcounts(image_orientation, Idof_edge);
            %             Idof_N       = circshift(Idof_N, 22); % shift 22 degree for fitting
            %************************* Gaussian fitting ********************
            centers      = [0 45 90 135];
            sigmas       = [8 8 8 8];
            amplitudes   = [1 1 1 1];
            %             centers      = [45 90 135];
            %             sigmas       = [8 8 8];
            %             amplitudes   = [1 1 1];
            numGaussians = 4;
            [amp, parameter, yhat] = fx_multi_Gaussianfitting(...
                numGaussians, Idof_edge(1:180), Idof_N, centers, sigmas, amplitudes);
            [C, amp_I]      = max(amp'.*parameter(2:2:end));
            m_fiber_angle   = parameter(amp_I*2-1);
            %             std_fiber_angle = amp(amp_I);
            %
            std_fiber_angle = std(image_orientation, 0, 'all','omitnan');
            disp(['mean:' num2str(m_fiber_angle) '. std:' num2str(std_fiber_angle)]);
            % ********** plot 1D distribution
            % normalize
            Idof_N = Idof_N / sum(Idof_N);
            % log
            figure('Name', ['angle_distribution_RT' '_' imagename]);
            set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
            h1 = bar(-89:90, Idof_N);
            xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
            %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
            set(gca, 'linewidth', 2);
            %         set(gca, 'YScale', 'log')
            %         legend([h1 h2], 'Angle distribution', 'Multiple Gaussian fitting', 'Interply track');
            %         ytickformat('percentage');
        end
        
        % ******************* RT ************************
        
        function obj = extract_local_orientation_RT_3D_zaxis(obj, PropertyName, radius, theta, z)
            % extract the local orientation image in the plane parallel to
            % surfaces
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % radius, theta: the ratio for the multi-resolution RT
            % z: the index along the z axis
            % ***
            % assign the info for saving
            obj.radius_RT             = radius;
            obj.theta_RT              = theta;
            obj.PropertyName_RT       = PropertyName;
            img_temp = abs(obj.(PropertyName));
            % start
            [lxc, lyc, ~]             = size(img_temp);
            Inplane_direction         = NaN(lxc, lyc, size(img_temp, 3));
            Inplane_direction_overall = NaN(size(obj.(PropertyName)));
            r_max                     = max(radius);
            %             ID_sum_all = NaN(length(distances_to_front), length(radius), length(theta));
            for i = 1:(length(z)-1)
                lower_index_bound = z(i+1);
                upper_index_bound = z(i);
                C_scan_inam_para_denoise = mean(img_temp(:, :, upper_index_bound:lower_index_bound-1), 3);
                % Radon transform
                %                 toc;
                for ii = 1+r_max: lxc - r_max
                    for jj = 1+r_max: lyc - r_max
                        for r_idx = 1: length(radius)
                            r = radius(r_idx);
                            center = [jj, ii];
                            % creat a round mask parameters
                            [~, anguler_1D, ~, ~] = fx_Radonto1Dangular_correct(C_scan_inam_para_denoise, center, r, theta);
                            %                             ID(ii, jj, r_idx, :) = anguler_1D;
                            % search for the max
                            [maxval, I] = max(anguler_1D);
                            % comment out to check the significance of the peak  normally needed
                            if maxval <= mean(anguler_1D) + std(anguler_1D) % there is or is not significant fibrous content.
                                continue;
                            end
                            Inplane_direction(ii, jj, upper_index_bound: lower_index_bound-1) = theta(I);
                        end
                    end
                end
                % end loop in one Cscan
                clc;
                fprintf('RT 3D progress: %0.2f%%\n',100*i/(length(z)-1));
            end
            obj.Inplane_direction_3D         = Inplane_direction;
            obj.Inplane_direction_3D_overall = Inplane_direction_overall;
            %             obj.ID_sum_overall               = ID_sum_all;
        end
        
        % ******************* ID **************************
        
        function obj = extract_local_orientation_3D_zaxis_allwl(obj, PropertyName, wavelength, orientation, z, K)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % z : index along z axis
            % K: the smoothing applied to the Gabor magnitude responses.
            % ***
            obj.wavelength_LG         = wavelength;
            obj.theta_LG              = orientation;
            obj.PropertyName_LG       = PropertyName;
            img_temp = abs(obj.(PropertyName));
            % start
            gaborArray                = gabor(wavelength, orientation, ...
                'SpatialFrequencyBandwidth', 1, 'SpatialAspectRatio', 0.7);
            row_ga                    = length(wavelength);
            col_ga                    = length(orientation);
            wl_max                    = round(max(wavelength)/2);
            Inplane_direction         = single(NaN(size(obj.(PropertyName))));
%             Inplane_direction_overall = single(NaN(size(obj.(PropertyName))));
%             ID_sum_all = NaN(length(z)-1, length(wavelength), length(orientation));
            %
            disp('gaborMagnitude calculation...');
            for i = 1:(length(z)-1)
                lower_index_bound = z(i+1);
                upper_index_bound = z(i);
                C_scan_inam_para_denoise = mean(img_temp(:,:,upper_index_bound:lower_index_bound-1), 3);
                [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_para_denoise, gaborArray);
                % normalized the magnitude.
                for k = 1:length(gaborArray)
                    BW         = gaborArray(k).SpatialFrequencyBandwidth;
                    sigmaX     = gaborArray(k).Wavelength / pi * sqrt(log(2)/2) * (2^BW+1) / (2^BW-1);
                    sigmaY     = sigmaX ./ gaborArray(k).SpatialAspectRatio;
                    % normalized the magnitude.
                    gabormag_k = gaborMagnitude(:, :, k);
                    gabormag_k = gabormag_k / (2 * sigmaX * sigmaY * pi);
                    % Gaussian filter for the GaborMag
                    if K==0 % choose to applie the filter
                        gaborMagnitude(:, :, k)  = gabormag_k;
                    else
                        sigma                    = 0.5*gaborArray(k).Wavelength;
                        invalid_part             = isnan(gabormag_k);
                        gabormag_k(invalid_part) = mean(gabormag_k, 'all', 'omitnan');
                        gabormag_k               = imgaussfilt(gabormag_k, K*sigma);
                        gabormag_k(invalid_part) = NaN;
                        gaborMagnitude(:, :, k)  = gabormag_k;
                    end
                end
                % search for the max
%                 ID_sum   = zeros(length(wavelength), length(orientation));
                [lx, ly] = size(gaborMagnitude(:, :, 1));
                for ii = wl_max + 1: lx - wl_max
                    for jj = wl_max + 1: ly - wl_max
                        ID     = squeeze(gaborMagnitude(ii, jj, :));
                        ID     = reshape(ID, [row_ga, col_ga]);
                        [C, I] = max(ID(:));
%                         if C <= mean(ID,'all') + std(ID, 0, 'all') % no significate maximum, thus skip
%                             continue;
%                         end
                        [~, indAng] = ind2sub(size(ID),I);
                        Inplane_direction(ii, jj, upper_index_bound: ...
                            lower_index_bound-1) = orientation(indAng);
%                         % sum up the ID
%                         ID_sum = ID_sum + ID;
                    end
                end
%                 % save the summary ID
%                 ID_sum_all(i, :, :) = ID_sum;
%                 [~, I]              = max(ID_sum(:));
%                 [~, indAng]         = ind2sub(size(ID_sum),I);
%                 Inplane_direction_overall(:, :, upper_index_bound: ...
%                     lower_index_bound) = orientation(indAng);
                % end loop in one Cscan
                clc;
                fprintf('ID 3D progress: %0.2f%%\n',100*i/(length(z)-1));
            end
            obj.Inplane_direction_3D_ID         = Inplane_direction;
%             obj.Inplane_direction_3D_overall_ID = Inplane_direction_overall;
%             obj.ID_sum_overall_ID               = ID_sum_all;
        end
        
        % ******************** 2d fft **********************
        
        function check2dfft(obj, z, PropertyName)
            % apply 2d fft to the C scan by z index.
            % depth.
            % z: z index. unit: data points
            % PropertyName: property to sue
            % ***
            % fillna
            C_scan_inam = abs(obj.(PropertyName)(:, :, z));
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            %
            % now make 2D fft of original image
            nfftx         = 2^nextpow2(m);
            nffty         = 2^nextpow2(n);
            % cannot solve the problem when nfftx!=nffty !
            nfftx         = max(nfftx, nffty);
            nffty         = max(nfftx, nffty);
            %
            fft2D         = fft2(C_scan_inam, nfftx, nffty);
            fft2D_shifted = fftshift(fft2D);
            % Gaussian Filter Response Calculation
            X              = -nfftx/2+1: nfftx/2;
            Y              = -nffty/2+1:nffty/2;
            % **** plot original
            cf = figure('Name', ['C_scan_amp' , '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 2);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            imagesc(ax, Yd, Xd, 20*log10(abs(C_scan_inam)/nfftx/nffty));
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            ax     = subplot(2, 1, 1);
            imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, 20*log10(abs(fft2D_shifted)/nfftx/nffty)); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
        end
        
        function obj = angulardistribution_2dfft(obj, zrange, PropertyName)
            % analyze angular distribution (2dfft) by z index (depth).
            % zrange: z index. unit: data points
            % PropertyName: property to sue
            % ***
            % fillna
            img_temp = abs(obj.(PropertyName));
            theta    = 0:180;
            angular_distribution_3D = nan(length(zrange), length(theta));
            for z_index = 1:length(zrange)
                z = zrange(z_index);
                C_scan_inam = img_temp(:,:,z);
                [m, n] = size(C_scan_inam);
                %
                % now make 2D fft of original image
                nfftx         = 2^nextpow2(m);
                nffty         = 2^nextpow2(n);
                % cannot solve the problem when nfftx!=nffty !
                nfftx         = max(nfftx, nffty);
                nffty         = max(nfftx, nffty);
                %
                fft2D         = fft2(C_scan_inam, nfftx, nffty);
                fft2D_shifted = fftshift(fft2D);
                % Gaussian Filter Response Calculation
                X              = -nfftx/2+1: nfftx/2;
                Y              = -nffty/2+1:nffty/2;     
                % ******  angular distribution
                theta    = 0:180;
                radii    = floor(min(nfftx, nffty)/2);
                [xx, yy] = meshgrid(-radii+1:radii, -radii+1:radii);
                mask     = (hypot(xx, yy)<=radii);
                fft2D_shifted = fft2D_shifted .* mask;
                %  *****　linear gaussian filter
                alpha = 30;
                w = gausswin(nfftx, alpha);
                w_2d =  repmat(w, 1, nffty);
                %　rotate
                angle = 135;
                w_2d_rot = imrotate(w_2d, angle, 'bilinear','crop');
                % fliter
                fft2D_shifted = fft2D_shifted.*(1-w_2d_rot);
                % **************
                [R,~] = radon(abs(fft2D_shifted), theta);
                % *************
                %                 % Display the transform.
                %                 imshow(R,[],'Xdata',theta,'Ydata',xp,'InitialMagnification','fit')
                %                 xlabel('\theta (degrees)')
                %                 ylabel('x''')
                %                 colormap(gca, jet), colorbar;
                % **********
                %                 % sum up
                angular_1D = sum(R(round(end/2-5):round(end/2+5), :), 1);
                % %                 "Gaussian" sum up
                %                 alpha = 30;
                %                 w = gausswin(size(R,1),alpha);
                %                 Rw = R.*w;
                %                 angular_1D = sum(Rw, 1);
                angular_1D = circshift(angular_1D, 22); % shift to remove edge effect
                angular_distribution_3D(z_index, :) = angular_1D/max(angular_1D); % by mean or max
                clc;
                fprintf('defect searching progress: %0.2f%%\n',100*z/max(zrange));
            end
            % **** plot distribution 3D ***
            cf = figure('Name', ['angular_distribution_3D']);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(1, 1, 1);
            imagesc(ax, (theta-22)-90, zrange, angular_distribution_3D);
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Normalized amp. (arb.)');
            xlabel('\fontname {times new roman} \theta (degree) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} z (point)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            obj.angular_distribution_fringerprint = angular_distribution_3D;
        end
        
        % ******************** 2D image-wise analysis *******************
        
        function obj = amplitude_rise_oneimage(obj, PropertyName, zrange, drop, filtertype)
            % determine thoverall size of image damagas
            % PropertyName: the name of the property in the object
            % zrange: the z axis range to find the max amplitudes
            % drop: the ratio of drop, -3 dB, -6 dB or .....
            % filtertype: type of the filter.
            img_temp = abs(obj.(PropertyName));
            img_temp = img_temp(:,:,zrange(1):zrange(2)); % select the range
            amp_max  = sum(img_temp, 3, 'omitnan');
            clear img_temp;
            % plot the 2d spectrum
            [m, n] = size(amp_max);
            % now make 2D fft of original image
            nfftx = 2^nextpow2(m);
            nffty = 2^nextpow2(n);
            nfftx = max(nfftx, nffty); % cannot tackle the issue of different nfftx and nffty.
            nffty = max(nfftx, nffty);
            %
            fft2D         = fft2(amp_max, nfftx, nffty);
            fft2D_shifted = fftshift(fft2D);
            X             = -nfftx/2+1: nfftx/2;
            Y             = -nffty/2+1:nffty/2;
            cf = figure('Name', ['2Dspectrum', '_', num2str(zrange)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(1, 1, 1);
            imagesc(X, Y, 20*log10(abs(fft2D_shifted)/m/n)); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %  % *** 2d fft filter
            switch filtertype{1}
                case char('LP')
                    prompt = 'What is the radius of the LP filter? ';
                    D_0 = input(prompt);
                    amp_max_filter = fx_lowpass_2dfft(amp_max, D_0, filtertype{2});
                case char('BP')
                    % Butterworth bandpass filter
                    prompt = 'What is the radius of the BP filter? ';
                    D_0 = input(prompt);
                    prompt = 'What is the width of the BP filter? ';
                    W = input(prompt);
                    if strcmp(filtertype{2}, 'Butterworth')
                        prompt = 'What is the order of the Butterworth filter? ';
                    end
                    n_bwf = input(prompt);
                    amp_max_filter = fx_bandpass_2dfft(amp_max, D_0, W, filtertype{2}, n_bwf);
                otherwise % no filter
                    amp_max_filter = amp_max;
            end
            % find the average peak amplitudes
            amp_max_mean = mean(amp_max_filter, 'all');
            log_amp_max  = 20*log10(amp_max_filter);
            groudtruth   = log_amp_max >= 20*log10(amp_max_mean) + drop;
            % display
            [m, n] = size(amp_max_filter);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            cf = figure('Name', ['defect_amplitduedrop', '_', num2str(drop)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 1);
            imagesc(ax, Yd, Xd, amp_max_filter);
            axis image;
            colormap(ax, jet);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % ground truth
            ax     = subplot(2, 1, 2);
            imagesc(ax, Yd, Xd, groudtruth);
            axis image;
            colormap(ax, gray);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Truth');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
        end
        
        function obj = amplitude_rise_3D(obj, PropertyName, zrange, drop, filtertype)
            % determine thoverall size of image damagas
            % PropertyName: the name of the property in the object
            % zrange: the z axis range to find the max amplitudes
            % drop: the ratio of drop, -3 dB, -6 dB or .....
            % filtertype: type of the filter.
            img_temp = abs(obj.(PropertyName));
            defect_groundtruth = boolean(zeros(size(img_temp)));
            %  % *** 2d fft filter
            switch filtertype{1}
                case char('LP')
                    prompt = 'What is the radius of the LP filter? ';
                    D_0 = input(prompt);
                case char('BP')
                    % Butterworth bandpass filter
                    prompt = 'What is the radius of the BP filter? ';
                    D_0 = input(prompt);
                    prompt = 'What is the width of the BP filter? ';
                    W = input(prompt);
                    if strcmp(filtertype{2}, 'Butterworth')
                        prompt = 'What is the order of the Butterworth filter? ';
                    end
                    n_bwf = input(prompt);
            end
            for z = zrange(1):zrange(2)
                amp_max = img_temp(:,:,z);
                % fill nan with zero for 2dfft2
                [row, col] = find(isnan(amp_max));
                amp_max(row, col) = mean(amp_max, 'all', 'omitnan');
                % filter
                switch filtertype{1}
                    case char('LP')
                        amp_max_filter = fx_lowpass_2dfft(amp_max, D_0, filtertype{2});
                    case char('BP')
                        amp_max_filter = fx_bandpass_2dfft(amp_max, D_0, W, filtertype{2}, n_bwf);
                    otherwise % no filter
                        amp_max_filter = amp_max;
                end
                % find the average peak amplitudes
                amp_max_filter(row, col) = NaN; % recover nan values
                amp_max_mean = mean(amp_max_filter, 'all', 'omitnan');
                log_amp_max  = 20*log10(amp_max_filter);
                groudtruth   = log_amp_max >= 20*log10(amp_max_mean) + drop; 
                defect_groundtruth(:,:,z) = groudtruth;
                % shadow effect
                [row,col] = find(groudtruth);
                img_temp(row, col, z+1:end) = NaN;                  
                clc;
                fprintf('defect searching progress: %0.2f%%\n',100*z/(zrange(2)-zrange(1)));
            end
            obj.defect_gt = defect_groundtruth;
        end
        
        function obj = amplitude_rise_oneimage_spatialfilter(obj, PropertyName, zrange, drop, sigma)
            % determine thoverall size of image damagas
            % PropertyName: the name of the property in the object
            % zrange: the z axis range to find the max amplitudes
            % drop: rise of amp in dB
            % sigma: sigma for Gaussian filter
            img_temp = obj.(PropertyName);
            inam_ex  = abs(img_temp(:,:,zrange(1):zrange(2))); % select the range    
            % ****
            inph_ex  = angle(img_temp);            
            infq_ex  = diff(unwrap(inph_ex, [], 3), 1, 3) /2 / pi * obj.fs;
            inph_ex  = sum(inph_ex(:,:,zrange(1):zrange(2)), 3, 'omitnan'); % select the range    
            infq_ex  = sum(infq_ex(:,:,zrange(1):zrange(2)), 3, 'omitnan'); % select the range
            amp_max  = sum(inam_ex, 3, 'omitnan');
            %             inph_ex  = angle(obj.(PropertyName));
            %            inph_temp = inph_ex(:,:,zrange(1):zrange(2)); % select the range
            clear img_temp;
            % Gaussian or median filter
            %             inph_temp = imgaussfilt(inph_temp, sigma);
            %             amp_max_filter = imgaussfilt(amp_max, sigma);
            %             amp_max_filter = medfilt2(amp_max, [sigma sigma]);
            
            %　2d wavelet filter
            amp_max = amp_max/max(amp_max, [], 'all')*127;
            [thr,sorh,keepapp] = ddencmp('den','wv', amp_max);
            amp_max_filter = wdencmp('gbl',amp_max,'sym4',2,thr,sorh,keepapp);
%             [amp_max_filter, ~, ~] = func_denoise_dw2d(amp_max/max(amp_max, [], 'all')*127);
            % find the average peak amplitudes
            amp_max_mean = mean(amp_max_filter, 'all', 'omitnan');
            log_amp_max  = 20*log10(amp_max_filter);
            groudtruth   = log_amp_max >= 20*log10(amp_max_mean) + drop;
            % imopen 
            se = strel('disk', 5);
            groudtruth = imopen(groudtruth, se);
            % display
            [m, n] = size(amp_max_filter);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            cf = figure('Name', ['defect_amplitduedrop', '_', num2str(drop)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 2, 1);
            imagesc(ax, Yd, Xd, amp_max_filter);
            axis image;
            colormap(ax, jet);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % inph
            ax     = subplot(2, 2, 2);
            imagesc(ax, Yd, Xd, inph_ex);
            axis image;
            colormap(ax, gray);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Truth');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % infq
            ax     = subplot(2, 2, 3);
            imagesc(ax, Yd, Xd, infq_ex);
            axis image;
            colormap(ax, gray);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Truth');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % ground truth
            ax     = subplot(2, 2, 4);
            imagesc(ax, Yd, Xd, groudtruth);
            axis image;
            colormap(ax, gray);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Truth');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
%             % edge detection 
%             BW1 = edge(amp_max_filter, 'sobel');
            
        end
        
        function obj = amplitude_rise_3D_spatialfilter(obj, PropertyName, zrange, drop, sigma)
            % determine thoverall size of image damagas
            % PropertyName: the name of the property in the object
            % zrange: the z axis range to find the max amplitudes
            % drop: rise of amp in dB
            % sigma: sigma for Gaussian filter
            img_temp = abs(obj.(PropertyName));
            %             inph_ex  = angle(obj.(PropertyName));
            defect_groundtruth = boolean(zeros(size(img_temp)));
            rear_I      = nan(size(img_temp, 1), size(img_temp, 2));
            rear_I_temp = rear_I;
            for z = zrange(1):zrange(2)
                amp_max = img_temp(:,:,z);
                % *** shadow effect
                %                 amp_max(~isnan(rear_I)) = mean(amp_max, 'all', 'omitnan');
                %                 amp_max(~isnan(rear_I)) = nan;
                % *** *
                %  Gaussian filter
                amp_max_filter = imgaussfilt(amp_max, sigma);
                %                 %　medfilt
                %                 amp_max_filter = medfilt2(amp_max, [sigma sigma]);
                %                 % meanfilter
                %                 h = 1/9*ones(9,1);
                %                 H = h*h';
                %                 amp_max_filter = filter2(H,amp_max);
%                 %　2d wavelet filter
%                 amp_max = amp_max/max(amp_max, [], 'all')*127;
%                 [thr,sorh,keepapp] = ddencmp('den','wv', amp_max);
%                 amp_max_filter = wdencmp('gbl',amp_max,'sym4',2,thr,sorh,keepapp);
                % find the average peak amplitudes
                amp_max_mean = mean(amp_max_filter, 'all', 'omitnan');
                log_amp_max  = 20*log10(amp_max_filter);
                if drop >=0
                    groundtruth   = log_amp_max >= 20*log10(amp_max_mean) + drop;
                else
                    groundtruth   = log_amp_max < 20*log10(amp_max_mean) + drop;
                end
                % imopen:　remove open gap
                se = strel('disk', 5);
                groundtruth = imopen(groundtruth, se);
                % imclose：remove strutrural gap
                groundtruth = imclose(groundtruth, se);
                %
                defect_groundtruth(:,:,z) = groundtruth;
                %                 if sum(groudtruth, 'all')>0
                rear_I_temp(groundtruth) = z;
                rear_I(rear_I~=0) = min(rear_I(rear_I~=0), rear_I_temp(rear_I~=0), 'omitnan') ; % rear_I ~=0
                rear_I(rear_I==0) = rear_I_temp(rear_I==0); % rear_I == 0
                %                 end
                clc;
                fprintf('defect searching progress: %0.2f%%\n',100*z/(zrange(2)-zrange(1)));
            end
            obj.defect_gt = defect_groundtruth;
            obj.rear_I = rear_I;
        end

        function obj = rear_filter(obj)
            rear_I_temp = obj. rear_I;
%             rear_I_temp(isnan(rear_I_temp)) = 650;
            [XDEN,~,~] = func_denoise_dw2d(rear_I_temp);
            obj. rear_I = XDEN;
        end
        
         function obj = search_rear_AS(obj, PropertyName, zrange)
            % determine thoverall size of image damagas
            % PropertyName: the name of the property in the object
            % zrange: the z axis range to find the max amplitudes
            % drop: the ratio of drop, -3 dB, -6 dB or .....
            % filtertype: type of the filter.
            img_temp = obj.(PropertyName);
            img_temp = img_temp(:,:,zrange(1):zrange(2)); % select the range  
            amp_max  = max(img_temp, [], 3, 'omitnan');
            temp_cp  = obj.c_p;
            temp_cp  = temp_cp(:,:,zrange(1):zrange(2)); % select the range  
            % 
            rear_I_temp = obj.rear_I;
            figure, imagesc(rear_I_temp);
            [m, n] = size(amp_max);
            % ************* stage 1: useing c_p and inph to correct the rear surface (no Generalization)***************** %
            inph_ex  = max(temp_cp, [], 3);
%             inph_ex = std(temp_cp, [], 3);
%             orthosliceViewer(temp_cp, 'Colormap', jet, 'ScaleFactors', [1 1 1]);
            %             close all;
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
%             figure, imagesc(Yd, Xd, inph_ex);
            groundtruth = inph_ex < mean(inph_ex, 'all', 'omitnan');
            % imclose：remove strutrural gap 
            se = strel('disk', 13);
            groundtruth = imclose(groundtruth, se);
%             figure, imagesc(groundtruth);
            % ************* phase analysis **********
            inph     = angle(img_temp);
            %             figure, plot(squeeze(inph(30, 30, :)));
            inph_cos        = cos(inph);
            inph_cos_filter = imgaussfilt3(inph_cos, 3, 'padding', 'replicate');
%             figure, plot(squeeze(inph_cos_filter(20, 20, :)));
%             rear_I_nondamage = rear_I_temp;
            for i = 1:m
                for j = 1:n
                    if ~groundtruth(i, j) % non-damage point
                        thershold = -0.3;
                        if find(inph_cos_filter(i ,j, :)<thershold, 1)
                            k = find(inph_cos_filter(i ,j, :)<thershold, 1); % returns the first n indices corresponding to the nonzero elements in X.
                            rear_I_temp(i ,j) = k + zrange(1);
                        end
                    end
                end
                clc;
                fprintf('%d / %d\n', i, m);
            end
            figure, imagesc(rear_I_temp);
            % ************* stage 2: to be filled****
            obj.rear_I = rear_I_temp;
         end
        
        %　** not yet corrected
        function obj = correct_rear_I_byinph(obj)
            % correct rear surface by instantanteous phase info.
            inph_ex     = angle(obj.img_hil);
            inph_unwarp = unwrap(inph_ex, [], 3);
            infq_ex     = diff(unwrap(inph_ex, [], 3), 1, 3) /2 / pi * obj.fs;
            temp_cp     = obj.c_p;
            %             orthosliceViewer(inph_ex);
            rear_I_temp  = obj.rear_I;
            rear_I_mean  = round(mean(rear_I_temp, 'all', 'omitnan'));
            front_I_temp = obj.front_I;
            front_I_mean = round(mean(front_I_temp, 'all', 'omitnan'));
            % index the 3d inph
            phase_rear        = single(nan(size(rear_I_temp)));
            infq_rear         = single(nan(size(rear_I_temp)));
            phase_front       = single(nan(size(rear_I_temp)));
            inph_unwarp_rear  = single(nan(size(rear_I_temp)));
            inph_unwarp_front = single(nan(size(rear_I_temp)));
            temp_cp_rear      = single(nan(size(rear_I_temp)));
            for i = 1:size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    % if front index is out of boundary, set mean
                    if isnan(front_I_temp(i, j)) || round(front_I_temp(i, j)) <= 0 || round(front_I_temp(i, j)) > size(inph_ex, 3)
                        phase_front(i, j) = inph_ex(i, j, front_I_mean);
                        inph_unwarp_front(i, j) = inph_unwarp(i, j, front_I_mean);
                    else
                        phase_front(i, j) = inph_ex(i, j, round(front_I_temp(i, j)));
                        inph_unwarp_front(i, j) = inph_unwarp(i, j, round(front_I_temp(i, j)));
                    end
                    % if rear index is out of boundary, set mean
                    if isnan(rear_I_temp(i, j)) || round(rear_I_temp(i, j)) > size(inph_ex, 3) || round(rear_I_temp(i, j)) <= 0
                        phase_rear(i, j)       = inph_ex(i, j, rear_I_mean);
                        infq_rear(i, j)        = infq_ex(i, j, rear_I_mean);
                        inph_unwarp_rear(i, j) = inph_unwarp(i, j, rear_I_mean);
                        temp_cp_rear(i, j)     = temp_cp(i, j, rear_I_mean);
                    else
                        phase_rear(i, j)  = inph_ex(i, j, round(rear_I_temp(i, j)));
                        infq_rear(i, j)  = infq_ex(i, j, round(rear_I_temp(i, j)));
                        inph_unwarp_rear(i, j)  = inph_unwarp(i, j, round(rear_I_temp(i, j)));
                        temp_cp_rear(i, j)     = temp_cp(i, j, round(rear_I_temp(i, j)));
                    end
                end
            end
%             figure,imagesc(temp_cp_rear);
%             figure,imagesc(rear_I_temp);
            
            orthosliceViewer(infq_ex);
            
            % correct by inph
            % **** the inph difference should be close to +-pi
            temp_inph_diff = mod((inph_unwarp_rear - inph_unwarp_front), pi);   
            figure,imagesc(temp_inph_diff);
            %             temp_inph_diff = phase_rear;
            %             temp_inph_diff(phase_rear>phase_front) =  ...
            %                 phase_rear(phase_rear>phase_front) - phase_front(phase_rear>phase_front);
            %             temp_inph_diff(phase_rear<=phase_front) =  ...
            %                 phase_front(phase_rear<=phase_front) - phase_rear(phase_rear<=phase_front);
            %             temp_inph_diff = mod(temp_inph_diff, pi);
            %             %
            %             temp = phase_rear<pi/2 & phase_rear>-pi/2;
            %             rear_I_temp(temp_inph_diff<pi/2) = nan;
            % **** the inph and infq should show continuty on the rear surface
            % local entrpyfilt;
            nhood = true(9);
%             J_inph      = entropyfilt(phase_rear, nhood);
%             figure, imagesc(phase_rear);
            J_infq = entropyfilt(infq_rear, nhood);   
%             J_infq_mask = J_infq < mean(J_infq, 'all');
%             figure, imagesc(J_infq_mask);
%             inph_mask = phase_rear>pi/1.5 | phase_rear < -pi/1.5;
%             figure, imagesc(temp_inph_diff);
            %
            %             rear_I_temp(J_inph>mean(J_inph, 'all')) = nan;
            rear_I_temp(J_infq>mean(J_infq, 'all')) = nan;
%             rear_I_temp(temp_inph_diff>pi*3/2 | temp_inph_diff<pi/2) = nan;
            % assign maximum depth to nan
            rear_I_temp(isnan(rear_I_temp)) = size(inph_ex, 3);
            obj.rear_I = rear_I_temp;
        end
        
        
        % ******************** cluster ***************
        
        function obj = knn_search_rear_I(obj, K, P, nbins)
            % search k nearest neighbors and display the hist of distance
            % K: number of the neighbors
            % parameters of the Minkowski distance
            % nbins: number of bins for the hist
            rear_I_temp      = obj.rear_I;
            [x, y]           = size(rear_I_temp);
            [x_mesh, y_mesh] = meshgrid(1:x, 1:y);
            %
            img_temp    = obj.img_hil;  
            inam_ex     = abs(img_temp); 
            inph_ex     = angle(obj.img_hil);
            inph_unwarp = unwrap(inph_ex, [], 3);
            infq_ex     = diff(unwrap(inph_ex, [], 3), 1, 3) /2 / pi * obj.fs;
            inph_unwarp_2d = nan(x, y);
            inam_2d = nan(x, y);
            infq_2d = nan(x, y);
            for i = 1:x
                for j = 1:y
                    % if front index is out of boundary, set mean
                    if isnan(rear_I_temp(i, j)) || round(rear_I_temp(i, j))>size(inph_unwarp, 3) || round(rear_I_temp(i, j))<=0
                        continue;
                    else
                        inph_unwarp_2d(i, j) = inph_unwarp(i, j, round(rear_I_temp(i, j)));
                        inam_2d(i, j) = inam_ex(i, j, round(rear_I_temp(i, j)));
                        infq_2d(i, j) = infq_ex(i, j, round(rear_I_temp(i, j)));
                    end
                end
            end
%             rear_scatter = [x_mesh(:) y_mesh(:) rear_I_temp(:)];
            rear_scatter = [rear_I_temp(:) inam_2d(:) inph_unwarp_2d(:) infq_2d(:)];
            rear_scatter_norm = normalize(rear_scatter, 2);
%             scatter3(rear_scatter(:, 1), rear_scatter(:, 2), rear_scatter(:, 3) ...
%                 , 5, rear_scatter(:, 3));
            
            [cIdx_temp, cD_temp] = knnsearch(rear_scatter_norm, rear_scatter_norm, 'K', K, 'NSMethod', 'kdtree', 'Distance', 'Minkowski', 'P', P);
            % hist the distances
            figure();
            histogram(mean(cD_temp, 2), nbins);
            set(gca,'YScale','log');
            xlabel('\fontname {times new roman} Average distance to k nearest neighbors (arb.) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Count', 'fontsize', 16);
            % discard the isolated points.
            mask = mean(cD_temp, 2)<=1e-7;
            figure,
            scatter3(x_mesh(mask), y_mesh(mask), rear_I_temp(mask) ...
                , 3, rear_I_temp(mask));
        end
    end
    
end

