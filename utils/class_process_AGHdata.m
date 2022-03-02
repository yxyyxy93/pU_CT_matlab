classdef class_process_AGHdata < class_process_RFdata
    % sub class of class_process_RFdata
    properties
    end
    
    methods
        function obj = class_process_AGHdata(filename)
            %
            obj = obj@class_process_RFdata(filename);
        end
        
        % ************************* preprocess *************************
        
        function obj = read_data_AGH(obj, path_load, n_files)
            % loading raw data
            % read one file to get the data size
            [B, scanSet.noA, scanSet.noP] = loadBin([path_load '0avg.bin']);
            
            % allocatin  g space
            img_origin = zeros(n_files, scanSet.noA, scanSet.noP);
            
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
            
            [lx, ly, lz] = size(img_origin);
            img_hilbert = zeros(lx, ly, lz);
            for i = 1:lx
                for j = 1:ly
                    img_hilbert(i, j, :) = hilbert(img_origin(i, j, :));
                end
                clc
                fprintf('Hilbert progress: %0.2f%%\n',100*i/lx)
            end
            obj.img = single(img_origin);
            obj.img_hil = single(img_hilbert);
        end
        
        function obj = normalize_timeAxis(obj)
            img_temp    = obj.img_hil;
            [lx, ly, ~] = size(img_temp);
            for i = 1:lx
                for j = 1:ly
                    img_temp(i, j, :) = img_temp(i, j, :)/max(abs(img_temp(i, j, :)));
                end
                clc;
                fprintf('normalization progress: %0.2f%%\n',100*i/lx);
            end%;
            obj.img_hil = img_temp;
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
            Inplane_direction         = NaN(size(obj.(PropertyName)));
            Inplane_direction_overall = NaN(size(obj.(PropertyName)));
            ID_sum_all = NaN(length(z)-1, length(wavelength), length(orientation));
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
                ID_sum   = zeros(length(wavelength), length(orientation));
                [lx, ly] = size(gaborMagnitude(:, :, 1));
                for ii = wl_max + 1: lx - wl_max
                    for jj = wl_max + 1: ly - wl_max
                        ID     = squeeze(gaborMagnitude(ii, jj, :));
                        ID     = reshape(ID, [row_ga, col_ga]);
                        [C, I] = max(ID(:));
                        if C <= mean(ID,'all') + std(ID, 0, 'all') % no significate maximum, thus skip
                            continue;
                        end
                        [~, indAng] = ind2sub(size(ID),I);
                        Inplane_direction(ii, jj, upper_index_bound: ...
                            lower_index_bound-1) = orientation(indAng);
                        % sum up the ID
                        ID_sum = ID_sum + ID;
                    end
                end
                % save the summary ID
                ID_sum_all(i, :, :) = ID_sum;
                [~, I]              = max(ID_sum(:));
                [~, indAng]         = ind2sub(size(ID_sum),I);
                Inplane_direction_overall(:, :, upper_index_bound: ...
                    lower_index_bound) = orientation(indAng);
                % end loop in one Cscan
                clc;
                fprintf('ID 3D progress: %0.2f%%\n',100*i/(length(z)-1));
            end
            obj.Inplane_direction_3D_ID         = Inplane_direction;
            obj.Inplane_direction_3D_overall_ID = Inplane_direction_overall;
            obj.ID_sum_overall_ID               = ID_sum_all;
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
    end
end

