classdef class_process_RFdata
    % Addme
    % process the cleaned RF raw volumeric data
    % a python script has been recommended to preprocess the data
    % then a cleaned data are produced including:
    % img_pre: preprocessed img
    % img_hil: analytical signal of preprocessed img
    % front: preprocessed front surface index
    % rear: preprocessed rear surface index
    
    properties
        data
        fns
        img
        img_hil
        front
        rear
        fs % Hz
        fx % / m
        fy % / m
        % surface
        front_I  % the new front index after process
        rear_I % the new rear index after process
        front_I_pre % without extra fillnan and smooth
        rear_I_pre % without extra fillnan and smooth
        rear_I_virtual % the vitual rear surface: ignoring the
        mask_wrong_front % a mask the same size as front_I to record the wrong front_I; 1: wrong
        % structural tensor
        c_p % the planar anisotropy metrics
        angle_x %  acos(eigenvector(1)) - pi/2;
        angle_y %  acos(eigenvector(2)) - pi/2;
        angle_z %  acos(eigenvector(3)) - pi/2;
        %
        img_hil_filter % the filtered AS (by logGabor filter)
        img_hil_noise % noise added by obj = addnoise(obj, snr)
        % ply track
        mask_plytrack
        inph_unwphase
        errors_pos
        est_thickness
        % the object of reference signal: 'class_reference_signal'
        reference
        refer_Ascan_ori
        refer_Ascan
        refer_Ascan_aligned
        % deconvolution
        img_WienerDeconv_B
        img_WienerDeconv_AR_B
        img_SparseDeconv_B
        img_WienerDeconv
        img_WienerDeconv_AR
        % properties for clustering
        layer_counts
        layer_counts_new
        X_scatter
        cIdx
        cD
        row
        col
        dep
        % in-plane
        C_scan_inam
        C_scan_inam_plywise
        C_scan_inam_denoise
        RT_ID
        RT_ID_gyy
        rear_mask
        Inplane_direction_3D
        Inplane_direction_3D_overall
        Inplane_direction_3D_overall_layer
        ID_sum_overall
        Inplane_direction_3D_modifyRT
        Inplane_direction_3D_overall_modifyRT
        Inplane_direction_3D_multireso_modifyRT
        ID_sum_overall_modifyRT
        % logGabor fitler
        gaborMag
        gaborPha
        frow
        fcol
        Inplane_direction_3D_ID
        Inplane_direction_3D_overall_ID
        ID_sum_overall_ID
        Inplane_direction_3D_ID_mean
        Inplane_direction_3D_ID_std
        % curvelet
        Curvelet_ID
        Curvelet_Coef
        % info of in-plane orientation extraction
        radius_RT
        theta_RT
        distances_to_front_RT
        PropertyName_RT
        wavelength_LG;
        theta_LG;
        distances_to_front_LG;
        PropertyName_LG;
        % one piexl info.
        ID_pixel
        % color define
        inst_amp_color   = [53 192 117]/255
        inst_phase_color = [0 172 239]/255
    end
    
    methods
        % ***************** definition, load data **************
        function obj = class_process_RFdata(filename)
            % filename: path and name of the mat. data including: img_pre, img_hil, front, rear
            if strcmp(filename(end-3:end), '.mat')
                obj.data = load(filename);
                obj.fns  = fieldnames(obj.data);
                disp(obj.fns);
            else
                disp('no data');
            end
        end
        
        function obj = read_data(obj)
            obj.img     = permute(obj.data.img_pre, [2 3 1]);
            obj.img_hil = permute(obj.data.img_hil, [2 3 1]);
            obj.front   = obj.data.front;
            obj.rear    = obj.data.rear;
            % new surfaces
            obj.front_I = obj.front;
            obj.rear_I = obj.rear;
        end
        
        function obj = read_origin_data(obj)
            img_origin = permute(obj.data.img_pre, [2 3 1]);
            [lx, ly, lz] = size(img_origin);
            img_hilbert = zeros(lx, ly, lz);
            for i = 1:lx
                for j = 1:ly
                    img_hilbert(i, j, :) = hilbert(img_origin(i, j, :));
                end
                disp(i);
            end
            obj.img = single(img_origin);
            obj.img_hil = single(img_hilbert);
        end
        
        function obj = read_data_fromMatpreprocess(obj)
            % read origin data from the preprocessed dataset by Matlab
            img_temp = obj.data.oimg_uint8;
            img_hil_temp = NaN(size(img_temp));
            for ii = 1: size(img_temp, 1)
                for jj = 1:size(img_temp, 2)
                    img_hil_temp(ii, jj, :) = hilbert(img_temp(ii, jj, :));
                end
            end
            obj.img = img_temp;
            obj.img_hil = img_hil_temp;
        end
        
        function obj = loadSettings(obj, settinig_file)
            % read the settings from the .xlsx file with the same name as
            % .mat: fs, fx, fy...
            % the .xlsx file is produced by preprocess_Cscan_AS\main.py
            [Samp, headertext] = xlsread(settinig_file, 'Settings');
            obj.fs          = str2double(headertext{9, 1}); % read the samping frequency
            if isnan(obj.fs)
                obj.fs = Samp;
            end
            obj.fx          = 1e3 / str2double(headertext{23, 1}); % read the x step and turn it into fx: mm -> 1 / m
            obj.fy          = 1e3 / str2double(headertext{25, 1}); % read the y step and turn it into fy: mm -> 1 / m
            % ! the sequcuse of x and y steps in the file are not confirmed
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
        
        function obj = shift_A_scan(obj, shift_len)
            % shift each A scan in one z diretion
            % shift_len: length of data points to shift
            img_temp     = obj.img;
            img_hil_temp = obj.img_hil;
            front_I_temp = obj.front_I;
            rear_I_temp  = obj.rear_I;
            % shift
            img_temp     = circshift(img_temp, shift_len ,3);
            img_hil_temp = circshift(img_hil_temp, shift_len ,3);
            front_I_temp = front_I_temp + shift_len;
            rear_I_temp  = rear_I_temp + shift_len;
            % save
            obj.img      = img_temp;
            obj.img_hil  = img_hil_temp;
            obj.front_I  = front_I_temp;
            obj.rear_I   = rear_I_temp;
        end
        
        function obj = time_varying_gain(obj, ratio)
            % apply a time_varying gain
            img_temp     = obj.img;
            img_hil_temp = obj.img_hil;
            [m, n, z]    = size(img_temp);
            for i = 1:m
                for j = 1:n
                    img_temp(i, j, :)     = squeeze(img_temp(i, j, :))...
                        .* linspace(1, z*ratio, z);
                    img_hil_temp(i, j, :) = squeeze(img_hil_temp(i, j, :))...
                        .* linspace(1, z*ratio, z);
                end
            end
            obj.img     = img_temp;
            obj.img_hil = img_hil_temp;
        end
        % ******************** display ****************
        function show_hilbert_Ascan(obj, x, y)
            % display the decomposed analytical signal of A scan
            % x, y: index of the A scan
            % reuse the fx_showAS
            ascan = squeeze(obj.img(x, y, :));
            ascan_hil = hilbert(ascan);
            fx_showAS((1:length(ascan)) / obj.fs, ascan_hil, obj.fs);
        end
        
        function show_Ascan_inam_peaks(obj, x, y, MinPD, MinPH, PropertyName)
            % demonstrate the inam of a ascan
            % find the peaks in the inam
            % x: x index
            % y: y index
            % MinPD: MinPeakDistance for findpeaks
            % MinPH: MinPeakHeight for findpeaks
            % PropertyName: the dataset to use
            % ***
            % reset front surfaces
            obj.front_I = obj.front;
            % calculate the inam
            % The sequence is changed unexpectedly. Transverse the matrix
            inam   = abs(obj.(PropertyName));
            inph   = angle(obj.(PropertyName));
            %
            t      = (1:size(inam, 3));
            A_inam = squeeze(inam(x, y, :));
            A_inph = squeeze(inph(x, y, :));
            figure;
            ca     = subplot(2, 1, 1);
            plot(t, A_inam, 'linewidth', 2);
            %             xlabel('\fontname {times new roman} Time(\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Inst. amp.', 'fontsize', 16);
            hold on;
            % find peaks
            [pks, locs, ~, ~] = findpeaks(A_inam, t, 'MinPeakDistance', MinPD, 'MinPeakHeight', MinPH);
            % select 2 maximum peaks
            [~, I]            = sort(pks);
            locs_max          = locs(I);
            plot(t(locs_max), A_inam(locs_max), 'rv');
            hold on;
            ca                = subplot(2, 1, 2);
            plot(t, A_inph, 'linewidth', 2);
            %             xlabel('\fontname {times new roman} Time(\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Inst. amp.', 'fontsize', 16);
            hold on;
            plot(t(locs_max), A_inph(locs_max), 'rv');
            hold on;
        end
        
        function obj = show_Cscan(obj, z, PropertyName)
            % show the C scan by z index in the depth
            % z: z index. unit: data points
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % In the image, the unit is transfered to us.
            %*********** origin
            %             figure('Name', ['C_scan_origin' , '_', num2str(z)]);
            %             C_scan = obj.(PropertyName)(:, :, z);
            %             ax = subplot(1, 1, 1);
            %             X = (1: size(C_scan, 1)) / obj.fx * 1e3;
            %             Y = (1: size(C_scan, 2)) / obj.fy * 1e3;
            %             imagesc(ax, X,Y, C_scan / obj.fx * 1e6)
            %             hold on;
            %             h = colorbar;
            %             %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            %             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            %             set(ax, 'fontsize', 16);
            %             set(ax, 'linewidth', 1.5);
            %             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(gca, 'linewidth', 2);
            %
            figure('Name', ['C_scan_amp' , '_', num2str(z)]);
            C_scan = abs(obj.(PropertyName)(:, :, z));
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan, 2)) / obj.fy * 1e3;
            C_scan = medfilt2(C_scan, [10 10]);
            imagesc(ax, Y, X, C_scan)
            hold on;
            h = colorbar; colormap(jet);
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
%             %
%             figure('Name', ['C_scan_inph' , '_', num2str(z)]);
%             C_scan = angle(obj.(PropertyName)(:, :, z));
%             ax = subplot(1, 1, 1);
%             X = (1: size(C_scan, 1)) / obj.fx * 1e3;
%             Y = (1: size(C_scan, 2)) / obj.fy * 1e3;
%             imagesc(ax, Y, X, C_scan / obj.fx * 1e6)
%             hold on;
%             h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
%             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
%             set(ax, 'fontsize', 16);
%             set(ax, 'linewidth', 1.5);
%             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
%             ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
%             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
%             set(gca, 'linewidth', 2);
            % save the Cscan
            obj.C_scan_inam  = C_scan;
        end
        
        function obj = addnoise(obj, snr)
            % add noise on the origin signal
            % snr: snr of the added noise
            [lx, ly, lz] = size(obj.img_hil);
            img_noiseadded = zeros(lx, ly, lz);
            for i = 1:lx
                for j = 1:ly
                    img_noiseadded(i, j, :) =  awgn(obj.img_hil(i, j, :), snr, 'measured');
                end
            end
            obj.img_hil_noise =  img_noiseadded;
        end
        
        function [Ascan, t_space, fss]= demo_Ascan(obj, x, y, fig_subname)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % fig_subname: the subname for the fig
            Ascan = squeeze(obj.img_hil_filter(x, y, :));
            % remove the direct component
            Ascan   = Ascan - mean(Ascan);
            fss     = obj.fs;
            t_space = (1:length(Ascan)) / fss * 1e6;
            % plot the inam and the record signal
            cf = figure('Name', strcat('Ascan_x', num2str(x), '_y', num2str(y), fig_subname) );
            set(cf, 'Position', [0, 0, 1000, 600], 'color', 'white');
            ca = subplot(1, 1, 1);
            plot(ca, t_space, real(Ascan), 'linewidth', 2);
            hold on;
            plot(ca, t_space, abs(hilbert(Ascan)), 'linewidth', 2);
            hold;
            legend({'Recorded signal', 'Envelope of the signal'});
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            %             ca = subplot(2, 1, 2);
            %             plot(ca, (1:length(Ascan)) / obj.fs, Ascan, 'linewidth', 2);
            %             ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            %             xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            %             legend({'Recorded', 'Inst. Amp.'});
            %             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(gca, 'linewidth', 2);
            %             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(gca, 'linewidth', 2);
            %             ylim([min(Ascan)/ 5, max(Ascan)/ 5]);
        end
        
        function demo_Bscan_interply(obj, B_type, index, Bwin, fig_subname)
            % demonstrate the B_scan of the img
            % show inph_ex
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % Bwin: the indexes of the first and the last points in Bscan
            % image
            % fig_subname: the subname for the fig
            inph_visual = obj.img;
            x = (1: Bwin(end)-Bwin(1)+1)  / obj.fx * 1e3;
            y = (1: Bwin(end)-Bwin(1)+1)  / obj.fy * 1e3;
            z = (1: size(inph_visual, 3)) / obj.fs * 1e6;
            cf = figure('Name', ['Bscan' '_', B_type, '_', num2str(index), fig_subname]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            %             ax = subplot(1, 1, 1);
            if (B_type == 'x')
                % B scan image
                B_scan = squeeze(inph_visual(Bwin, index, :));
            elseif (B_type == 'y')
                % B scan image
                B_scan = squeeze(inph_visual(index, Bwin, :));
            end
            % choose the 'x' axis
            if (B_type == 'y')
                x = y;
            end
            % remove the direct component
            B_scan = B_scan - mean(B_scan, 2);
            imagesc(x, z, abs(hilbert(B_scan')));
            hold on;
            % figure setup
            colormap(gray);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            cl = caxis;
            caxis(cl / 3);
        end
        
        function show_inaminph_3D(obj, xslice, yslice, zslice)
            % demonstrate the 3d original results: inam, inph
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % angle_compens: the compensation for angle
            % ***
            inph_ex  = angle(obj.img_hil);
            inam_ex = abs(obj.img_hil);
            infq_ex   = diff(unwrap(inph_ex, [], 3), 1, 3) /2 / pi * obj.fs;
            if isempty(inph_ex)
                inph_ex = obj.Inplane_direction_3D_ID;
            end
            y = (0: size(inam_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inam_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inam_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
            % ******************** show inam ***********
            cf        = figure('Name', ['3d_inam_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax        = subplot(1, 1, 1);
            [X, Y, Z] = meshgrid(x, y, z);
            h         = slice(ax, X, Y, Z, inam_ex , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap gray;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            xlim([x(1) x(end)]);
            ylim([y(1) y(end)]);
%             % *********** surfaecs xslice **** !uncomment for displaying the surfaces
%             x_idx = xslice * obj.fx / 1e3;
%             y_idx = yslice * obj.fy / 1e3;
%             h1    = scatter3(xslice * ones(1, size(inph_ex, 1)), ...
%                 y, obj.front_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
%                 2, 'red', 'filled', 'DisplayName','Front surface');
%             hold on;
%             h2    = scatter3(xslice * ones(1, size(inph_ex, 1)), ...
%                 y, obj.rear_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
%                 2, 'magenta', 'filled', 'DisplayName','Rear surface');
%             hold on;
%             % select the yslice
%             scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
%                 obj.front_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
%                 2, 'red', 'HandleVisibility','off');
%             hold on;
%             scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
%                 obj.rear_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
%                 2, 'magenta', 'HandleVisibility','off');
%             hold on;
%             legend([h1 h2], 'Front surface', 'Rear surface');
%             zlim([1 5.5])
            view([15 65 40]);
            % ******************** show inph ***********
            cf = figure('Name', ['3d_inph_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
            [X, Y, Z] = meshgrid(x, y, z);
            h         = slice(ax, X, Y, Z, inph_ex , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap gray;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Phase (rad.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            xlim([x(1) x(end)]);
            ylim([y(1) y(end)]);
%             % surfaecs xslice v**** !uncomment for displaying the surfaces
%             x_idx = xslice * obj.fx / 1e3;
%             y_idx = yslice * obj.fy / 1e3;
%             h1 = scatter3(xslice * ones(1, size(inph_ex, 1)), ...
%                 y, obj.front_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
%                 2, 'red', 'filled', ...
%                 'DisplayName','Front surface');
%             hold on;
%             h2 = scatter3(xslice * ones(1, size(inph_ex, 1)), ...
%                 y, obj.rear_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
%                 2, 'magenta', 'filled', ...
%                 'DisplayName','Rear surface');
%             hold on;
%             % surfaecs yslice
%             scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
%                 obj.front_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
%                 2, 'red', 'HandleVisibility','off');
%             hold on;
%             scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
%                 obj.rear_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
%                 2, 'magenta', 'HandleVisibility','off');
%             hold on;
%             legend([h1 h2], 'Front surface', 'Rear surface');
            view([15 65 40]);
            % ******************** show infq ***********
            cf = figure('Name', ['3d_infq_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            z         = (0: size(infq_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
            [X, Y, Z] = meshgrid(x, y, z);
            h         = slice(ax, X, Y, Z, infq_ex , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap gray;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Phase (rad.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            xlim([y(1) y(end)]);
            ylim([x(1) x(end)]);
            % surfaecs xslice
            x_idx = xslice * obj.fx / 1e3;
            y_idx = yslice * obj.fy / 1e3;
            h1 = scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                y, obj.front_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
                2, 'red', 'filled', ...
                'DisplayName','Front surface');
            hold on;
            h2 = scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                y, obj.rear_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
                2, 'magenta', 'filled', ...
                'DisplayName','Rear surface');
            hold on;
            % surfaecs yslice
            scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                obj.front_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
                2, 'red', 'HandleVisibility','off');
            hold on;
            scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                obj.rear_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
                2, 'magenta', 'HandleVisibility','off');
            hold on;
            legend([h1 h2], 'Front surface', 'Rear surface');
            view([15 65 40]);
        end
        
        function show_Raw_filtered_3D(obj, xslice, yslice, zslice)
            % demonstrate the 3d original results and the filtered results
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % angle_compens: the compensation for angle
            % ***
            inam_ex     = real(obj.img_hil);
            inam_filter = real(obj.img_hil_filter);
            %
            y = (0: size(inam_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inam_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inam_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
            % ******************** show raw ***********
            cf        = figure('Name', ['3d_raw_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax        = subplot(1, 1, 1);
            [X, Y, Z] = meshgrid(x, y, z);
            h         = slice(ax, X, Y, Z, inam_ex , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            xlim([x(1) x(end)]);
            ylim([y(1) y(end)]);
            view([15 65 40]);
            % ******************** show filter ***********
            cf = figure('Name', ['3d_filtered_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            z         = (0: size(inam_filter, 3) - 1) / obj.fs * 1e3 * 3000/2;
            [X, Y, Z] = meshgrid(x, y, z);
            h         = slice(ax, X, Y, Z, inam_filter , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            xlim([x(1) x(end)]);
            ylim([y(1) y(end)]);
            view([15 65 40]);
        end
        
        function C_scan_inam_temp = display_slice(obj, C_scan_index, PropertyName)
            % demonstrate the slice of certain content
            % C_scan_index: xy index of the slice xy 
            % PropertyName: name of the 3d volumetric data
            inam = obj.(PropertyName);
            if strcmp(PropertyName, 'angle_z')
                inam = (inam + pi/2) / pi * 180;
            elseif strcmp(PropertyName, 'Inplane_direction_3D_ID')
                angle_compens = 0;
                inam = mod(inam + angle_compens, 180) - 90;
            end
            C_scan_inam_temp  = NaN([size(inam, 1), size(inam, 2)]);
            for i = 1: size(inam, 1)
                for j = 1:size(inam, 2)
                    % take both front and back surfaces into consideration
                    index = C_scan_index(i, j);
                    C_scan_inam_temp(i, j) = inam(i, j, round(index));
                end
            end
%            filter the profile and asign the inam
            figure('Name', ['Cscan_slice_' , PropertyName, num2str(round(mean(C_scan_index, 'all')))]);
            ax = subplot(1, 1, 1);
            X = (1: size(inam, 2)) / obj.fx * 1e3;
            Y = (1: size(inam, 1)) / obj.fy * 1e3;
            h = imagesc(ax, X, Y, C_scan_inam_temp);
            set(h, 'AlphaData', 1-isnan(C_scan_inam_temp)) % NaN = white
            hold on;
            h = colorbar;
            colormap jet;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} ');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        % ***************** surface determination ****************
        function [A, A0, z] = statistic_front_backechoes(obj, window)
            % window: the window selected for computation
            % caculate the average amplitude of the front- and back-surface
            % echoes as well as the distance.
            % A: average amplitude of back-surface echo
            % A0: average amplitude of front-surface echo
            % z: average distance(data points) between the front- and back-surface echoes
            x1       = window(1);
            x2       = window(2);
            y1       = window(3);
            y2       = window(4);
            t1       = window(5);
            t2       = window(6);
            inam_img = abs(obj.img_hil(x1:min(end, x2), y1:min(end, y2), t1: min(end, t2)));
            front_I_temp = obj.front_I(x1:min(end, x2), y1: min(end, y2));
            rear_I_temp  = obj.rear_I(x1:min(end, x2), y1:min(end, y2));
            amp_front    = NaN(size(front_I_temp));
            amp_rear     = NaN(size(rear_I_temp));
            A_scan       = squeeze(mean(inam_img, [1 2]));
            for i = 1:size(inam_img, 1)
                for j = 1:size(inam_img, 2)
                    % if front index is out of boundary, set mean
                    if isnan(front_I_temp(i, j)) || round(front_I_temp(i, j)) <= 0 ...
                            || round(front_I_temp(i, j))>size(inam_img, 3)
                        continue;
                    else
                        amp_front(i, j)   = inam_img(i, j, round(front_I_temp(i,j)));
                    end
                    % if rear index is out of boundary, set mean
                    if isnan(rear_I_temp(i,j)) || round(rear_I_temp(i,j)) > size(inam_img,3) ...
                            || round(rear_I_temp(i,j)) <= 0
                        continue;
                    else
                        amp_rear(i, j)    = inam_img(i, j, round(rear_I_temp(i, j)));
                    end
                end
            end
            A0 = mean(amp_front, [1,2], 'omitnan');
            A  = mean(amp_rear, [1,2], 'omitnan');
            z  = mean(rear_I_temp-front_I_temp, [1,2], 'omitnan');
        end
        
        function obj = find_front_amp(obj,  MinPD, MinPH, PropertyName, max_len, threshold_delamination)
            % determining the front surface index by amp.
            % obj.rear_i is also intially determined here
            % The sequence is changed unexpectedly. Transverse the matrix
            % MinPD: MinPeakDistance for findpeaks
            % MinPH: MinPeakHeight for findpeaks
            % PropertyName: 'img_hil' or 'img_hil_filter' ... used for the inam.
            % max_len: the end point in the signal
            % threshold_delamination: the threshold to determine the delaminiation
            % ****
            % calculate the inam
            % The sequence is changed unexpectedly. Transverse the matrix
            inam         = abs(obj.(PropertyName));
            inph         = angle(obj.(PropertyName));
            % the time domain
            t            = (1:max_len);
            front_I_temp = NaN(size(inam, 1), size(inam, 2));
            rear_I_temp  = NaN(size(inam, 1), size(inam, 2));
            %
            for i = 1:size(inam, 1)
                for j = 1:size(inam, 2)
                    A_inam            = squeeze(inam(i, j, 1:max_len));
                    A_inph            = squeeze(inph(i, j, 1:max_len));
                    [pks, locs, ~, ~] = findpeaks(A_inam, t, 'MinPeakDistance', MinPD, 'MinPeakHeight', MinPH);
                    %                     % sort the peaks
                    %                     [~, I] = sort(pks);
                    % the front-wall echo is largest
                    if length(pks)>1
                        % the first echo should be front
                        front_I_temp(i, j) = min(locs);
%                         % derive the front_I by phase
%                         A_inph_locs        = abs(A_inph(locs) - A_inph(min(locs)));
%                         [~, I_phase]       = sort(A_inph_locs);
                        % find the fisrt element larger than the front echo
                        pkf_thres = pks(1) * threshold_delamination;
                        if length(find(pks>pkf_thres, 2)) > 1
                            locs_pfk_thres    = locs(find(pks>pkf_thres, 2));
                            rear_I_temp(i, j) = locs_pfk_thres(2);
%                         elseif locs(I_phase(1)) < locs(I_phase(end))
%                             rear_I_temp(i, j) = locs(end);
                            % front_I_temp(i, j) = locs(I_phase(1));
                        else
%                             [~, I]            = sort(pks, 'descend');
%                             rear_I_temp(i, j) = locs(I(2));
                            rear_I_temp(i, j) = locs(end);
                        end
                    elseif length(pks)==1
                        front_I_temp(i, j)    = locs(1);
                    else
                        front_I_temp(i, j)    = NaN;
                    end
                end
                disp([num2str(i) '/' num2str(size(inam, 1))]);
            end
            front_I_temp    = fx_inpaint_nans(front_I_temp, 0);
            obj.front_I     = front_I_temp;
            obj.rear_I      = rear_I_temp;
            obj.front_I_pre = front_I_temp;
            obj.rear_I_pre  = rear_I_temp;
        end
        
        function obj = find_front_amp_alpha(obj,  MinPD, MinPH, PropertyName, max_len, alpha, A_ratio)
            % determining the front surface index by amp.
            % obj.rear_i is also intially determined here
            % define the delamination by the attenuated threshold!
            % The sequence is changed unexpectedly. Transverse the matrix
            % MinPD: MinPeakDistance for findpeaks
            % MinPH: MinPeakHeight for findpeaks
            % PropertyName: 'img_hil' or 'img_hil_filter' ... used for the inam.
            % max_len: the end point in the signal
            % threshold_delamination: the threshold to determine the delaminiation
            % Ab: the amplitude of the back-wall echo
            % ****
            % calculate the inam
            % The sequence is changed unexpectedly. Transverse the matrix
            global min_pks % set the minimus peak of the first echo
            inam         = abs(obj.(PropertyName));
%             inph         = angle(obj.(PropertyName));
            % the time domain
            t            = (1:max_len);
            front_I_temp = NaN(size(inam, 1), size(inam, 2));
            rear_I_temp  = NaN(size(inam, 1), size(inam, 2));
            %
            for i = 1:size(inam, 1)
                for j = 1:size(inam, 2)
                    A_inam            = squeeze(inam(i, j, 1:max_len));
                    [pks, locs, ~, ~] = findpeaks(A_inam, t, 'MinPeakDistance', MinPD, 'MinPeakHeight', MinPH);
                    %                     % sort the peaks
                    %                     [~, I] = sort(pks);
                    % the front-wall echo is largest
                    if length(pks)>1
                        % the first echo should be front
                        % front_I_temp(i, j) = min(locs);    
                        % set a threshold for front-surface echo as well
                        if find(pks>min_pks, 1, 'first')
                            loc_findex         = find(pks>min_pks, 1, 'first');
                            front_I_temp(i, j) = locs(loc_findex);
                        else 
                            continue;
                        end
                        % remove the front-wall echo
                        pks  = pks(loc_findex+1:end);
                        locs = locs(loc_findex+1:end);
                        if isempty(locs)
                            continue;
                        end
                        Af   = A_inam(front_I_temp(i, j));
                        % find the delamination by the attenuated threshold
                        pks_cmp = pks' - Af * A_ratio * exp(- alpha*(locs - front_I_temp(i, j)));
                        if find(pks_cmp>0, 1)
                            rear_I_temp(i, j) = locs(find(pks_cmp>0, 1));
                        else
                            rear_I_temp(i, j) = locs(end);
                        end
                    elseif length(pks)==1
                        front_I_temp(i, j) = locs(1);
                    else
                        front_I_temp(i, j) = NaN;
                    end
                end
                clc;
                disp([num2str(i) '/' num2str(size(inam, 1))]);
            end
            front_I_temp    = fx_inpaint_nans(front_I_temp, 0);
            obj.front_I     = front_I_temp;
            obj.rear_I      = rear_I_temp;
            obj.front_I_pre = front_I_temp;
            obj.rear_I_pre  = rear_I_temp;
        end
        
        function obj = find_front_amp_alpha_Ascan(obj,  MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, x, y)
            % determining the front surface index by amp.
            % obj.rear_i is also intially determined here
            % define the delamination by the attenuated threshold!
            % The sequence is changed unexpectedly. Transverse the matrix
            % MinPD: MinPeakDistance for findpeaks
            % MinPH: MinPeakHeight for findpeaks
            % PropertyName: 'img_hil' or 'img_hil_filter' ... used for the inam.
            % max_len: the end point in the signal
            % alpha: compensation factor
            % x,y: the index for Ascan
            % ****
            % calculate the inam
            % The sequence is changed unexpectedly. Transverse the matrix
            ascan        = obj.(PropertyName);
            ascan        = squeeze(ascan(x, y, 1:max_len));
            % inph         = angle(obj.(PropertyName));
            % the time domain
            t            = (1:max_len);
            t_domain     = t / obj.fs*1e6;
            %
            A_inam            = abs(ascan);
            [pks, locs, ~, ~] = findpeaks(A_inam, t, 'MinPeakDistance', MinPD, 'MinPeakHeight', MinPH);
            %                     % sort the peaks
            %                     [~, I] = sort(pks);
            % the front-wall echo is largest
            % the first echo should be front
            % front_I_temp(i, j) = min(locs);
            % set a threshold for front-surface echo as well
            global min_pks
            min_pks = 0.2;
            if find(pks>min_pks, 1, 'first')
                loc_findex    = find(pks>min_pks, 1, 'first');
                front_I_temp  = locs(loc_findex);
            end
            % remove the front-wall echo
            pks  = pks(loc_findex+1:end);
            locs = locs(loc_findex+1:end);
            Af   = A_inam(front_I_temp);
            % find the delamination by the attenuated threshold
            pks_cmp = pks' - Af * A_ratio * exp(-alpha*(locs - front_I_temp));
            if find(pks_cmp>0, 1)
                rear_I_temp = locs(find(pks_cmp>0, 1));
            else
                rear_I_temp = locs(end);
            end
            Gate = Af * A_ratio * exp(-alpha*(t - front_I_temp));
            % plot in fig.
            cf = figure('Name', ['Ascan_', 'TimeGate_xy_', num2str(x), '_', num2str(y)]);
            set(cf, 'Position', [0, 0, 800, 300], 'color', 'white');
            ca = subplot(1, 1, 1);
            h1 = plot(t_domain, ascan, 'k-', 'linewidth', 2);
            hold on;
            h2 = plot(t_domain, A_inam, '-', 'linewidth', 2, 'Color', [1 0.5 0]); % orange color
            hold on;
            h3 = plot(t_domain(t>front_I_temp), Gate(t>front_I_temp), '--', 'linewidth', 2);
            hold on;
            plot(t_domain(front_I_temp), Af, 'ro', 'linewidth', 2);
            hold on;
            plot(t_domain(rear_I_temp), A_inam(rear_I_temp), 'ro', 'linewidth', 2);
            hold on;
            legend([h1 h2 h3], 'Original signal', 'Instantaneous amplitude', 'Dynamic amplitude gate');
            xlabel('\fontname {times new roman} Time (\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Amp. (arb.)', 'fontsize', 16);
            set(ca, 'fontsize', 16);  
            set(ca, 'fontname', 'Times new roman');
            set(ca, 'linewidth', 1.5);
            xlim([t_domain(1) t_domain(end)]);
        end
        
        function obj = smooth_rear_I(obj, ~)
            % fill and smooth the rear_I
            % modefilter applied!!
            obj.front_I_pre = obj.front_I;
            obj.rear_I_pre  = obj.rear_I;
            smoothing_size  = 3;
            front_I_temp    = fx_inpaint_nans(obj.front_I, 5);
            front_I_temp    = padarray(front_I_temp, [smoothing_size smoothing_size], 'replicate');
            front_I_temp    = round(colfilt(front_I_temp, ...
                [smoothing_size smoothing_size], 'sliding', @mode));
            obj.front_I     = front_I_temp(smoothing_size+1:end-smoothing_size, ...
                smoothing_size+1:end - smoothing_size);
%             delete the odd points
%             rear_I_mean   = mean(obj.rear_I, 'all', 'omitnan');
            rear_I_temp   = fx_inpaint_nans(obj.rear_I, 5);
            rear_I_temp   = padarray(rear_I_temp, [smoothing_size smoothing_size], 'replicate');
            rear_I_temp   = round(colfilt(rear_I_temp, ...
                [smoothing_size smoothing_size], 'sliding', @min)); % take the minimum
            obj.rear_I    = rear_I_temp(smoothing_size+1:end-smoothing_size, ...
                smoothing_size+1:end - smoothing_size);
%             obj.rear_I  = round(modefilt(obj.rear_I, [3, 3], 'symmetric'));
        end
        
        function obj = recover_surface(obj)
            % recover the origin surface without smooth
            obj.front_I = obj.front_I_pre;
            obj.rear_I = obj.rear_I_pre;
        end
        
        function obj = find_rear_amp(obj, error)
            % determining the rear surface, as well as the vitual rear
            % surface for the delamination part
            % obj.rear_I is already intially determined by 'obj.find_front_amp'
            % error: the interval of the index to determine the delamination part
            rear_I_nan  = obj.rear_I;
            % calculate the Most frequent values in array
            rear_I_mean = round(mode(obj.rear_I, 'all'));
            for i = 1:size(obj.front_I, 1)
                for j = 1:size(obj.front_I, 2)
                    % find the misleading front surface by
                    % 'obj.wrong_front_mask'
                    if obj.mask_wrong_front(i, j) == 1
                        obj.rear_I(i, j) = obj.front_I(i, j);
                    end
                    % determine the delamination part
                    if abs(obj.rear_I(i, j) - rear_I_mean) > error
                        rear_I_nan(i, j) = NaN;
                    end
                end
            end
            rear_I_nan = round(medfilt2(rear_I_nan, [4, 4], 'symmetric'));
            % interpolate the rear_I_nan for the virtual rear_I
            rear_I_virtual_unfilter = round(fx_inpaint_nans(rear_I_nan, 5));
            obj.rear_I_virtual = round(medfilt2(rear_I_virtual_unfilter, [3, 3], 'symmetric'));
        end
        
        function show_surfaces(obj, filename)
            % display the surfaces: phase, amp., index
            % calculate the inam and inph
            % The sequence is changed unexpectedly. Transverse the matrix
            % filename: the path to save the fig
            inph = atan2(imag(obj.img_hil), real(obj.img_hil));
            %             inph = permute(inph, [2 3 1]);
            inam = abs(obj.img_hil);
            %             inam = permute(inam, [2 3 1]);
            % phase calculation
            phase_front  = zeros(size(obj.front_I));
            amp_front    = zeros(size(obj.front_I));
            phase_rear   = zeros(size(obj.front_I));
            amp_rear     = zeros(size(obj.front_I));
            %             phase_front_sum = 0;
            front_I_mean = round(mean(obj.front_I, 'all', 'omitnan'));
            rear_I_mean  = round(mean(obj.rear_I, 'all', 'omitnan'));
            %             rear_I_mean = round(mean(mean(obj.rear_I)));
            for i = 1:size(inam, 1)
                for j = 1:size(inam, 2)
                    % if front index is out of boundary, set mean
                    if isnan(obj.front_I(i, j)) || round(obj.front_I(i, j)) <= 0 || round(obj.front_I(i, j)) > size(inph, 3)
                        phase_front(i, j) = inph(i, j, front_I_mean);
                        amp_front(i, j)   = inam(i, j, front_I_mean);
                    else
                        phase_front(i, j) = inph(i, j, round(obj.front_I(i, j)));
                        amp_front(i, j)   = inam(i, j, round(obj.front_I(i, j)));
                    end
                    % if rear index is out of boundary, set mean
                    if isnan(obj.rear_I(i, j)) || round(obj.rear_I(i, j)) > size(inph, 3) || round(obj.rear_I(i, j)) <= 0
                        phase_rear(i, j)  = inph(i, j, rear_I_mean);
                        amp_rear(i, j)    = inam(i, j, rear_I_mean);
                    else
                        phase_rear(i, j)  = inph(i, j, round(obj.rear_I(i, j)));
                        amp_rear(i, j)    = inam(i, j, round(obj.rear_I(i, j)));
                    end
                end
            end
            %             phase_rear_mean = phase_front_mean + pi; % rad
            % visual the phase_front amp_front
            cf = figure('name', filename);
            ca = subplot(3, 1, 1);
            imagesc(ca, phase_front); colorbar;
            title('The phase of the surface');
            ca = subplot(3, 1, 2);
            imagesc(ca, amp_front); colorbar;
            title('The amp. of the surface');
            ca = subplot(3, 1, 3);
            imagesc(ca, obj.front_I); colorbar;
            xlabel('y');
            ylabel('x');
            title('The index. of the surface');
            % comment this if necessary
%             % save the fig
%             saveas(cf, filename, 'bmp');
%             saveas(cf, filename, 'fig');
            % rear surface
            figure, ca = subplot(3, 1, 1);
            imagesc(ca, phase_rear); colorbar;
            title('The phase of the surface');
            ca = subplot(3, 1, 2);
            imagesc(ca, amp_rear); colorbar;
            title('The amp. of the surface');
            ca = subplot(3, 1, 3);
            imagesc(ca, obj.rear_I); colorbar;
            title('The index. of the surface');
            xlabel('y');
            ylabel('x');
        end
        
        function show_rear_surface(obj)
            % display the rear and virtual rear surfaces: index
            figure, ca = subplot(2, 1, 1);
            imagesc(ca, obj.rear_I); colorbar;
            title('The index. of the rear surface');
            xlabel('y');
            ylabel('x');
            ca = subplot(2, 1, 2);
            imagesc(ca, obj.rear_I_virtual); colorbar;
            title('The index. of the virtual rear surface');
            xlabel('y');
            ylabel('x');
        end
        
        function show_inph_3d(obj, PropertyName, xslice, yslice, zslice)
            % demonstrate the 3d instantaneous phase
            % show inph_ex
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            %                         xslice =  100 / obj.fx * 1e3;
            %             yslice = 50 / obj.fy * 1e3;
            %             zslice = [];
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            inph = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            inph_ex = fx_ExcludeData(inph, 'f', obj.front, 'r', obj.rear);
            y = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            inph_visual = inph_ex;
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            ax = subplot(1, 1, 1);
            h = slice(ax, X, Y, Z, inph_visual, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap gray;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. phase(rad.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            set(gca,'zDir','reverse'); % get time domain upside down
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
        end
        
        function obj = cut_edges(obj, window, x_step, y_step)
            % cut the edges of the dataset, and display the processed
            % figure
            % window: [x1, x2, y1, y2, t1, t2]
            % x_step, y_step: the steps to change the scanning intervel
            x1 = window(1);
            x2 = window(2);
            y1 = window(3);
            y2 = window(4);
            t1 = window(5);
            t2 = window(6);
            obj.img     = obj.img(x1: x_step: min(end, x2), y1: y_step: min(end, y2), t1: min(end, t2));
            obj.img_hil = obj.img_hil(x1: x_step: min(end, x2), y1: y_step: min(end, y2), t1: min(end, t2));
            obj.front   = obj.front(x1: x_step: min(end, x2), y1: y_step: min(end, y2)) - t1;
            obj.rear    = obj.rear(x1: x_step: min(end, x2), y1: y_step: min(end, y2)) - t1;
            obj.front_I = obj.front_I(x1: x_step: min(end, x2), y1: y_step: min(end, y2)) - t1;
            obj.rear_I  = obj.rear_I(x1: x_step: min(end, x2), y1: y_step: min(end, y2)) - t1;
            
            % change index step
            obj.fx      = obj.fx / y_step; % fx is related to y_step
            obj.fy      = obj.fy / x_step;
        end
        
        function obj = perfect_rear_surface(obj)
            % make the rear_surface to the end of the 3D volumetric dataset
            obj.rear_I = size(obj.img, 3) * ones(size(obj.img, 1), size(obj.img, 2));
            obj.rear_mask = [];
        end
        
        function obj = show_front_position(obj, winx, winy, rot_angle)
            % display the front surface profie
            % xslice, yslice: the index to slice the profile
            % winx, winy: the window indexes to display the data
            % rot_angle: angle for rotating the image
            front_profile = obj.front_I / obj.fs * 1e6 * 1500 / 2;
            % this is a rotation to fit the optical data!!
            front_profile = imrotate(front_profile, rot_angle, 'crop');          
            front_profile = front_profile(winx, winy);
            % detrend
            %Make the 2D data as 1D vector
            [xl, yl] = size(front_profile);
            [Y, X]   = meshgrid(1:yl, 1:xl);
            Xcolv    = X(:); % Make X a column vector
            Ycolv    = Y(:); % Make Y a column vector
            Zcolv    = front_profile(:); % Make Z a column vector
            Const    = ones(size(Xcolv)); % Vector of ones for constant term
            % find the coeffcients of the best plane fit
            Coefficients = [Xcolv Ycolv Const]\Zcolv; % Find the coefficients
            XCoeff = Coefficients(1); % X coefficient
            YCoeff = Coefficients(2); % X coefficient
            CCoeff = Coefficients(3); % constant term
            Z_p           = XCoeff * X + YCoeff * Y + CCoeff;
            front_profile = front_profile - Z_p;
            front_profile = -(front_profile - max(front_profile(:)));
            [xl, yl]      = size(front_profile);
            Z_f           = front_profile;
            % 2D profile image
            cf = figure('Name', 'surface_profile_top_view');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            imagesc((1:yl)/obj.fx*1e3, (1:xl)/obj.fy*1e3, Z_f);
            shading flat;
            % set(gca, 'ydir', 'reverse');
            hold on;
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Depth (\mum)');
            xlabel('X displacement (mm)');
            ylabel('Y displacement (mm)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            set(ax, 'fontname', 'Times new roman');
%             xticks(2:2:18);
            pbaspect([length(winy), length(winx), 1]);
            % surface lines plot  
            % x and y surface lines
            x_line = Z_f(round(end/2), :);
            y_line = Z_f(:, round(end/2));
            cf = figure('Name', 'surface_lines_xy');
            set(cf, 'Position', [0, 0, 1200, 600], 'color', 'white');
            % x
            ax     = subplot(2, 1, 1);
            x_axis = (1:yl)/obj.fy*1e3;
            plot(x_axis, x_line, 'linewidth', 2);
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Depth (mm)', 'fontsize', 16);
            % y
            ax = subplot(2, 1, 2);
            y_axis = (1:xl)/obj.fx*1e3;
            plot(y_axis, y_line, 'linewidth', 2);
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} Y displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Depth (mm)', 'fontsize', 16);
            save('profile_line.mat', 'x_axis', 'y_axis', 'x_line', 'y_line');
        end
        
        function obj = show_front_position_circle(obj, centers, radii)
            % display the front surface profie in a masked circle
            % xslice, yslice: the index to slice the profile
            % winx, winy: the window indexes to display the data
            front_profile = obj.front_I / obj.fs * 1e6 * 1500 / 2;
            % define rectangular window
            x1            = centers(2) - radii;
            x2            = centers(2) + radii;
            y1            = centers(1) - radii;
            y2            = centers(1) + radii;
            %
            [xl, yl]      = size(front_profile);
            [Y, X]        = meshgrid(1:yl, 1:xl);
            % detrend
            %Make the 2D data as 1D vector
            Xcolv  = X(:); % Make X a column vector
            Ycolv  = Y(:); % Make Y a column vector
            Zcolv  = front_profile(:); % Make Z a column vector
            Const  = ones(size(Xcolv)); % Vector of ones for constant term
            % find the coeffcients of the best plane fit
            Coefficients = [Xcolv Ycolv Const]\Zcolv; % Find the coefficients
            XCoeff = Coefficients(1); % X coefficient
            YCoeff = Coefficients(2); % X coefficient
            CCoeff = Coefficients(3); % constant term
            Z_p           = XCoeff * X + YCoeff * Y + CCoeff;
            front_profile = front_profile - Z_p;
            front_profile = -(front_profile - max(front_profile(:)));
%             Z_f           = front_profile;
            % crop a circle
            % Crop and asign the image
            inam_C_scan      = front_profile;
            inam_C_scan_mask = inam_C_scan(x1: x2, y1: y2);
            [xx, yy]         = meshgrid(-radii:radii, -radii:radii);
            mask             = (hypot(xx, yy)<=radii);
            inam_C_scan_mask = double(inam_C_scan_mask) .* mask;
            Z_f              = inam_C_scan_mask;
            % 2D profile image
            cf = figure('Name', 'surface_profile_top_view');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            imagesc((1:size(Z_f, 2))/obj.fx*1e3, (1:size(Z_f, 1))/obj.fy*1e3, Z_f);
            shading flat;
            % set(gca, 'ydir', 'reverse');
            hold on;
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Depth (\mum)');
            xlabel('X displacement (mm)');
            ylabel('Y displacement (mm)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
        end
        
        % ************** internal damage features imaging
        function obj = damage_imaging(obj)
            % plot the internal damge features in 3D
            back_surface = (obj.rear_I - obj.front_I) / obj.fs * 1e3 * 3000 / 2; 
            inam         = abs(obj.img_hil);
            [Y, X] = ndgrid(1:size(back_surface,1), 1:size(back_surface,2));
            Y            = Y / obj.fy * 1e3;
            X            = X / obj.fy * 1e3;
            amp_rear     = zeros(size(obj.front_I));
            % phase_front_sum = 0;
            rear_I_index = obj.rear_I;
            rear_I_mean  = round(mean(rear_I_index, 'all', 'omitnan'));
            % rear_I_mean = round(mean(mean(obj.rear_I)));
            for i = 1:size(inam, 1)
                for j = 1:size(inam, 2)
                    % if rear index is out of boundary, set mean
                    if isnan(obj.rear_I(i, j)) || round(obj.rear_I(i, j)) > size(inam, 3) || round(obj.rear_I(i, j)) <= 0
                        amp_rear(i, j)    = inam(i, j, rear_I_mean);
                    else
                        amp_rear(i, j)    = inam(i, j, round(rear_I_index(i,j)));
                    end
                end
            end
            % compute the amp. along depth
%             figure;
%             scatter(rear_I_index(:), amp_rear(:), 1, 'filled');
            % top view of amplitude
            cf = figure('Name', 'Internal_damage_features_amp_top_view');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            pcolor(X, Y, amp_rear);
            shading flat;
            hold on;
            colormap jet;
            h  = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            % top view of depth
            cf = figure('Name', 'Internal_damage_features_dep_top_view');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            pcolor(X, Y, back_surface);
            shading flat;
            % set(gca, 'ydir', 'reverse');
            hold on;
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Depth (mm)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            % angled view 
            cf = figure('Name', 'Internal_damage_features');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            scatter3(X(:), Y(:), back_surface(:), 3, amp_rear(:), 'filled');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Depth (mm)', 'fontsize', 16);
            set(gca,'zDir','reverse'); % get time domain upside down
            xlim([1 max(X(:))]);
            ylim([1 max(Y(:))]);
            zlim([0 6]);
            view([15 65 40]);
        end
        
        % ******************** filter ***********
        function obj = Filter_lowpass(obj, bandpassFreq, PropertyName)
            % Filter the 3d dataset by low-pass filter
            % bandpassFreq: passband freq. of the filter
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            [lx, ly, ~] = size(obj.(PropertyName));
            % save in the obj.img_hil_filter
            obj.img_hil_filter = zeros(size(obj.(PropertyName)));
            % define a low-pass filter and filter the signals
            [~, lpFilt] = lowpass(squeeze(obj.(PropertyName)(1, 1, :)), bandpassFreq, obj.fs, ...
                'ImpulseResponse', 'fir', 'Steepness', 0.8, 'StopbandAttenuation', 60);
            delay       = mean(grpdelay(lpFilt));
            for x = 1: lx
                for y = 1:ly
                    Ascan    = squeeze(obj.(PropertyName)(x, y, :));
                    Ascan_as = filter(lpFilt, Ascan);
                    Ascan_as = circshift(Ascan_as, -delay);
                    obj.img_hil_filter(x, y, :) = Ascan_as;
                end
                clc;
                disp([num2str(x), '/' , num2str(lx)]);
            end
        end
        
        function show_logGabor_Scaleogram(obj, f0, sigma, x, y)
            % demonstrate the Scaleogram of inph of various scales of
            % logGabor filter
            % f0: the array of the center frequencies of the filter
            % sigma: affects the bandwidth of the filter.
            % x, y: the x, y index of the 3d dataset
            Ascan = squeeze(obj.img(x, y, :));
            Scaleogram = zeros( length(f0), length(Ascan) );
            for i = 1:length(f0)
                % calculate the AS in frequency domain
                Fo = fft(Ascan);
                % transverse Fo for dot multiply
                Fo = Fo';
                fx_shift = fftshift(Fo);
                % Fo(-w) = conj(Fo(w));
                n = length(Ascan);
                fshift = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
                Fa = (1 + sign(fshift)) .* fx_shift;
                % add the logGabor filter here
                filter =  fx_1dLogGabor_filter(fshift, f0(i), sigma);
                Fa_filter = Fa .* filter;
                % not filtered, for debug
                % ifft to obtain the AS
                %                 Fa_ishift = ifftshift(Fa);  % back to the original transform output
                %                 as = ifft(Fa_ishift);
                %                 figure, plot(abs(as));
                % ifft to obtain the AS
                Fa_ishift = ifftshift(Fa_filter);  % back to the original transform output
                as = ifft(Fa_ishift);
                % the direction of 'as' is reversed somehow, needs to be
                % flipped and conjugated!
                Scaleogram(i, :) = conj(flip(as));
            end
            figure, imagesc( (1:length(Ascan)) / obj.fs, f0, angle(Scaleogram));
            colorbar;
            hold on;
            line([obj.front_I(x, y) / obj.fs; obj.front_I(x, y)  / obj.fs], [f0(1), f0(end)], ...
                'linestyle', '--', 'LineWidth', 2, 'Color', 'r');
            hold on;
            line([obj.rear_I(x, y) / obj.fs; obj.rear_I(x, y)  / obj.fs], [f0(1), f0(end)], ...
                'linestyle', '--', 'LineWidth', 2, 'Color', 'r');
            ylabel({['\fontname {times new roman}\fontsize {16}' 'Center frequency of filter (Hz)']});
            xlabel('\fontname {times new roman} Time(s)', 'fontsize', 16);
        end
        
        function demo_logGabor_plytrack(obj, x, y, f0, sigma)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % n_interface: the number of the interfacec we want to track
            % f0: the array of the center frequencies of the filter
            % sigma: affects the bandwidth of the filter.
            Ascan = squeeze(obj.img(x, y, :));
            front_I_interp = obj.front_I(x, y);
            rear_I_interp = obj.rear_I(x, y);
            % define a logGabor filter and calcualte the analytical signal
            % calculate the AS in frequency domain
            Fo = fft(Ascan);
            % transverse Fo for dot multiply
            Fo = Fo';
            fx_shift = fftshift(Fo);
            n = length(Ascan);
            fshift = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            Fa = (1 + sign(fshift)) .* fx_shift;
            % add the logGabor filter here
            filter =  fx_1dLogGabor_filter(fshift, f0, sigma);
            Fa_filter = Fa .* filter;
            % ifft to obtain the AS
            Fa_ishift = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as = conj(flip(Ascan_as));
            % plot the inam and the record signal
            cf = figure('Name', strcat('plytrack_filter', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(2, 1, 1);
            plot(ca, (1:length(Ascan_as)) / obj.fs, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, (1:length(Ascan_as)) / obj.fs, abs(Ascan_as), 'linewidth', 2)
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Filtered', 'Inst. Amp.'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % track the - pi / 2 by interplation on the unwraped data
            A_inph        = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            A_inph_unw    = unwrap(A_inph);
            phase_seq     = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
            I_phase_track = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                phase_seq) + front_I_interp;
            % drop the Nan in 'I_phase_track'
            I_phase_track = I_phase_track(~isnan(I_phase_track));
            % plot the inph, and the ply track of - pi / 2
            ca = subplot(2, 1, 2);
            plot(ca, (1:length(Ascan_as)) / obj.fs, angle(Ascan_as), 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, [front_I_interp rear_I_interp] / obj.fs , angle(Ascan_as(round([front_I_interp rear_I_interp]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            hold on;
            plot(ca, I_phase_track / obj.fs , angle(Ascan_as(round(I_phase_track))) , 'gd', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Inst. phase', 'Derived surface', 'Derived inter-ply'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function demo_logGabor_plytrack_inph(obj, x, y, f0, sigma, threshold)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % n_interface: the number of the interfacec we want to track
            % f0: the array of the center frequencies of the filter
            % sigma: affects the bandwidth of the filter.
            % threshold: threshold is a ratio to find the peaks
            Ascan                  = squeeze(obj.img(x, y, :));
            % define a logGabor filter and calcualte the analytical signal
            % calculate the AS in frequency domain
            Fo                     = fft(Ascan);
            % transverse Fo for dot multiply
            Fo                     = Fo';
            fx_shift               = fftshift(Fo);
            n                      = length(Ascan);
            fshift                 = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            Fa                     = (1 + sign(fshift)) .* fx_shift;
            % add the logGabor filter here
            filter                 = fx_1dLogGabor_filter(fshift, f0, sigma);
            Fa_filter              = Fa .* filter;
            % ifft to obtain the AS
            Fa_ishift              = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as               = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as               = conj(flip(Ascan_as));
            % find the front and back
            inam_ascan             = abs(Ascan_as);
            MinPeakHeight          = threshold * max(inam_ascan);
            MinPeakDistance        = length(inam_ascan) / 3;
            [~, Idx_walls]         = findpeaks(inam_ascan, 1:length(inam_ascan), ...
                'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
            front_I_interp         = Idx_walls(1);
            rear_I_interp          = Idx_walls(end);
            % plot the inam and the record signal
            %             % track the - pi/2 by interplation on the unwraped data
            Idx_resin              = linspace(front_I_interp, rear_I_interp, 25);
            %             [lct, ~, ~, ~, A_inph] = fx_track_byphase(Idx_resin, real(Ascan_as), threshold);
            % track the -pi/2 by distance and clustering
            %             [lct, inam_ascan, Ascan_phase] = fx_track_bycluster(Ascan_as, threshold, 23, Idx_resin(2:end-1));
            %             [lct, inam_ascan, Ascan_phase, errors_pos] = fx_track_bydistance(Ascan_as, threshold, Idx_resin);
            % track the -pi/2 by phase
            [lct, pct, errors_pos, inam_ascan, Ascan_phase] = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);     
            %             % debug
            %             figure, plot(Ascan_phase);
            %             hold on;
            %             scatter(round(lct_C), Ascan_phase(round(lct_C)));
            % plot the inph, and the ply track of - pi / 2
            % plot the inam and the record signal
            cf             = figure('Name', strcat('plytrack_filter', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca             = subplot(2, 1, 1);
            t_space        = (1:length(Ascan_as)) / obj.fs * 1e6; % us
            plot(ca, t_space, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, t_space, inam_ascan, 'linewidth', 2)
            hold on;
            plot(ca, t_space([front_I_interp rear_I_interp]), ...
                inam_ascan([front_I_interp rear_I_interp]), 'ro', 'linewidth', 2);
            hold on;
            plot(ca, t_space(lct([1 end])), inam_ascan(lct([1 end])), 'ro', 'linewidth', 1);
            hold on;
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Recorded signal', 'Inst. Amp.', 'Determined positions of surfaces'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            ca = subplot(2, 1, 2);
            plot(ca, t_space, Ascan_phase, 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, t_space(lct([1 end])), Ascan_phase(lct([1 end])), 'ro', 'linewidth', 1);
            hold on;
            plot(ca, lct(2:end-1)/obj.fs*1e6, pct(2:end-1), 'ms', 'linewidth', 1);
            hold on;
            plot(ca, t_space(round(Idx_resin(2:end-1))), Ascan_phase(round(Idx_resin(2:end-1))), 'gd', 'linewidth', 1);
            hold on;
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Inst. phase', 'Determined positions of surfaces', 'Determined positions of interplies', 'Theoretical positions of interplies'}); 
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % for debug
            figure,
            plot(errors_pos);
        end
        
        function demo_logGabor_plytrack_inph_v2(obj, x, y, f0, sigma, threshold)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % n_interface: the number of the interfacec we want to track
            % f0: the array of the center frequencies of the filter
            % sigma: affects the bandwidth of the filter.
            % threshold: threshold is a ratio to find the peaks
            Ascan                  = squeeze(obj.img(x, y, :));
            % define a logGabor filter and calcualte the analytical signal
            % calculate the AS in frequency domain
            Fo                     = fft(Ascan);
            % transverse Fo for dot multiply
            Fo                     = Fo';
            fx_shift               = fftshift(Fo);
            n                      = length(Ascan);
            fshift                 = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            Fa                     = (1 + sign(fshift)) .* fx_shift;
            % add the logGabor filter here
            filter                 = fx_1dLogGabor_filter(fshift, f0, sigma);
            Fa_filter              = Fa .* filter;
            % ifft to obtain the AS
            Fa_ishift              = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as               = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as               = conj(flip(Ascan_as));
            % find the front and back
            inam_ascan             = abs(Ascan_as);
            MinPeakHeight          = threshold * max(inam_ascan);
            MinPeakDistance        = length(inam_ascan) / 2;
            [~, Idx_walls]         = findpeaks(inam_ascan, 1:length(inam_ascan), ...
                'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
            front_I_interp         = round(Idx_walls(1));
            rear_I_interp          = round(Idx_walls(end));
            % plot the inam and the record signal
            %             % track the - pi/2 by interplation on the unwraped data
            Idx_resin              = linspace(front_I_interp, rear_I_interp, 25);
            %             [lct, ~, ~, ~, A_inph] = fx_track_byphase(Idx_resin, real(Ascan_as), threshold);
            % track the -pi/2 by distance and clustering
            %             [lct, inam_ascan, Ascan_phase] = fx_track_bycluster(Ascan_as, threshold, 23, Idx_resin(2:end-1));
            %             [lct, inam_ascan, Ascan_phase, errors_pos] = fx_track_bydistance(Ascan_as, threshold, Idx_resin);
            % track the -pi/2 by phase
            [lct, ~, errors_pos, inam_ascan, Ascan_phase] = fx_track_byphase(Idx_resin, Ascan_as, threshold);     
            %             % debug
            %             figure, plot(Ascan_phase);
            %             hold on;
            %             scatter(round(lct_C), Ascan_phase(round(lct_C)));
            % plot the inph, and the ply track of - pi / 2
            % plot the inam and the record signal
            cf             = figure('Name', strcat('plytrack_filter_v2_', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca             = subplot(2, 1, 1);
            t_space        = (1:length(Ascan_as)) / obj.fs * 1e6; % us
            plot(ca, t_space, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, t_space, inam_ascan, 'linewidth', 2)
            hold on;
            plot(ca, t_space([front_I_interp rear_I_interp]), ...
                inam_ascan([front_I_interp rear_I_interp]), 'ro', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Recorded signal', 'Inst. Amp.', 'Determined positions of surfaces'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            ca = subplot(2, 1, 2);
            plot(ca, t_space, Ascan_phase, 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, t_space(lct([1 end])), Ascan_phase(lct([1 end])), 'ro', 'linewidth', 1);
            hold on;
            plot(ca, t_space(round(Idx_resin(2:end-1))), Ascan_phase(round(Idx_resin(2:end-1))), 'gd', 'linewidth', 1);
            hold on;
            yl = yline(ca, Ascan_phase(lct(1))-pi/2,'--','LineWidth',2);
            yl.Color = [.80 0 .40];
            text(t_space(end), Ascan_phase(lct(1))-pi/2+pi/5, '\phi_0-\pi/2', 'Color', [.80 0 .40],'FontSize',16);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Inst. phase', 'Determined positions of surfaces', 'Theoretical positions of interplies'}); 
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % for debug
            figure,
            plot(errors_pos);
        end
        
        function demo_lowpass_plytrack(obj, x, y, n_interface)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % n_interface: the number of the interfacec we want to track
            % method: "cutoff"
            Ascan = squeeze(obj.img(x, y, :));
            front_I_interp = obj.front_I(x, y);
            rear_I_interp = obj.rear_I(x, y);
            % define a logGabor filter and calcualte the analytical signal
            % calculate the AS in frequency domain
            Fo = fft(Ascan);
            % transverse Fo for dot multiply
            Fo = Fo';
            fx_shift = fftshift(Fo);
            n = length(Ascan);
            fshift = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            Fa = (1 + sign(fshift)) .* fx_shift;
            % add the logGabor filter here
            filter = [ ];
            Fa_filter = Fa .* filter;
            % ifft to obtain the AS
            Fa_ishift = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as  = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as  = conj(flip(Ascan_as));
            % plot the inam and the record signal
            cf = figure('Name', strcat('plytrack_lowpassfilter', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(2, 1, 1);
            plot(ca, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, abs(Ascan_as), 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amp. (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Recorded', 'Inst. Amp.'})
            % track the - pi / 2 by interplation on the unwraped data
            A_inph     = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            A_inph_unw = unwrap(A_inph);
            phase_seq = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : ...
                (2 * n_interface - 2) * pi +  A_inph_unw(1) + 3 * pi / 2;
            I_phase_track = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                phase_seq) + front_I_interp;
            % plot the inph, and the ply track of - pi / 2
            ca = subplot(2, 1, 2);
            plot(ca, angle(Ascan_as), 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, I_phase_track, angle(Ascan_as(round(I_phase_track))) , 'bv', 'linewidth', 2);
            hold on;
            plot(ca, [front_I_interp rear_I_interp], angle(Ascan_as(round([front_I_interp rear_I_interp]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} ', 'fontsize', 16);
        end
        
        function demo_logGabor_plytrack_addnoise(obj, x, y, f0, sigma, snr_value)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % n_interface: the number of the interfacec we want to track
            % f0: the array of the center frequencies of the filter
            % sigma: affects the bandwidth of the filter.
            % snr: the snr of the added noise
            Ascan          = squeeze(obj.img(x, y, :));
            front_I_interp = obj.front_I(x, y);
            rear_I_interp  = obj.rear_I(x, y);
            % add noise
            Ascan          = awgn(Ascan, snr_value, 'measured');
            % define a logGabor filter and calcualte the analytical signal
            % calculate the AS in frequency domain
            Fo             = fft(Ascan);
            % transverse Fo for dot multiply
            Fo             = Fo';
            fx_shift       = fftshift(Fo);
            n              = length(Ascan);
            fshift         = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            Fa             = (1 + sign(fshift)) .* fx_shift;
            % add the logGabor filter here
            filter         = fx_1dLogGabor_filter(fshift, f0, sigma);
            Fa_filter      = Fa .* filter;
            % ifft to obtain the AS
            Fa_ishift      = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as       = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as       = conj(flip(Ascan_as));
            % calcualte the snr
            %             SNR = snr(real(Ascan_as), obj.fs);
            pn             = bandpower(real(Ascan_as), obj.fs, [10e6, obj.fs / 2.1]);
            ps             = bandpower(real(Ascan_as), obj.fs, [1e6, 10e6]);
            SNR            = 10*log10(ps/pn);
            disp("SNR:");
            disp(SNR);
            % % track the - pi / 2 by interplation on the unwraped data
            % A_inph        = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            % A_inph_unw    = unwrap(A_inph);
            % phase_seq     = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
            % I_phase_track = interp1(A_inph_unw, 1:length(A_inph_unw), ...
            %   phase_seq) + front_I_interp;
            % track the - pi / 2 by interplation on the wraped data
            Idx_resin           = linspace(front_I_interp, rear_I_interp, 25);
            threshold           = 0.05;
            [lct, pct, ~, ~, ~] = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);
            % drop the Nan in 'I_phase_track'
            I_phase_track       = lct(2:end-1);
            I_phase_track       = I_phase_track(~isnan(I_phase_track));
            P_phase_track       = pct(2:end-1);
            P_phase_track       = P_phase_track(~isnan(P_phase_track));
            % plot the inam and the record signal
            cf             = figure('Name', strcat('plytrack_filter', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca             = subplot(2, 1, 1);
            plot(ca, (1:length(Ascan_as))/obj.fs*1e6, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, (1:length(Ascan_as))/obj.fs*1e6, abs(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, [front_I_interp lct(end)]/obj.fs*1e6, abs(Ascan_as(round([front_I_interp lct(end)]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Filtered', 'Inst. Amp.'});
            legend('Location', 'best');
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % plot the inph, and the ply track of - pi / 2
            ca = subplot(2, 1, 2);
            plot(ca, (1:length(Ascan_as))/obj.fs*1e6, angle(Ascan_as), 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            % use extracted location of rear surface here.
            plot(ca, [front_I_interp lct(end)]/obj.fs*1e6, angle(Ascan_as(round([front_I_interp lct(end)]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            hold on;
            plot(ca, I_phase_track/obj.fs*1e6, P_phase_track, 'gd', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Inst. phase', 'Derived surface', 'Derived inter-ply'});
            legend('Location', 'northwest');
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        % ***************** structural tensor to extract angle ****************
        function obj = structural_tensor(obj, sigma1, sigma2, PropertyName, ds_rate)
            % apply structural tensor on the cleaned inph dataset
            % explanation about the angles:
            % anglex, y, z represent the angles with the axis directional vector of the planar like structure.
            % extract the angles of the the plane-structrual
            % sigma1: smoothing scale
            % sigma2: integration scale
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % ds_rate: the rate of downsampling
            % ***
            % time record
            tic;
            disp('Start structural tensor method');
            % The sequence is changed unexpectedly. Transverse the matrix
            inph = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            % downsampling
            n         = ds_rate;
            inph      = inph(1:n:end, 1:n:end, :);
            % 3D Gaussian filter
            inph_cos        = cos(inph);
            inph_sin        = sin(inph);
            inph_cos_filter = imgaussfilt3(inph_cos, sigma1, 'padding', 'replicate');
            inph_sin_filter = imgaussfilt3(inph_sin, sigma1, 'padding', 'replicate');
            % calculation of gradients
            [Gx_cos, Gy_cos, Gz_cos] = fx_imgradientxyz(inph_cos_filter, 'intermediate');
            [Gx_sin, Gy_sin, Gz_sin] = fx_imgradientxyz(inph_sin_filter, 'intermediate');
            Gx = inph_cos.*Gx_sin - inph_sin.*Gx_cos;
            Gy = inph_cos.*Gy_sin - inph_sin.*Gy_cos;
            Gz = inph_cos.*Gz_sin - inph_sin.*Gz_cos;
            % scale phase gradients in units of phaseper- millimetre to ensure correct angular measurements
            Gx = Gx * obj.fx;
            Gy = Gy * obj.fy;
            Gz = Gz * obj.fs * 3.7e-6 / 5.5e-3;
            % release the memory
            clear inph_cos;
            clear inph_sin;
            clear inph;
            clear inam;
            clear infq;
            clear inph_cos_filter;
            clear inph_sin_filter;
            clear Gx_cos;
            clear Gy_cos;
            clear Gz_cos;
            clear Gx_sin;
            clear Gy_sin;
            clear Gz_sin;
            % for debug
            % [~, ~] = fx_showBS(2, inph, num_bs, step_bs, data.front, data.rear, B_type, 'jet');
            %             [~, ~] = fx_showBS(3, Gx, num_bs, step_bs, data.front, data.rear, B_type, 'jet');
            %             [~, ~] = fx_showBS(4, Gy, num_bs, step_bs, data.front, data.rear, B_type, 'jet');
            %             [fig, ax] = fx_showBS(5, Gz, num_bs, step_bs, data.front, data.rear, B_type, 'jet');
            % form the structure-tensor
            ST    = fx_formStructureTensor(Gx, Gy, Gz);
            % downsample
            ST_in = zeros(size(ST));
            for i=1:3
                for j=1:3
                    ST_in(:, :, :, i, j) = imgaussfilt3(ST(:, :, :, i, j), sigma2);
                end
            end
            % time record
            timeElapsed = toc;
            disp(['form_tensor: ', num2str(timeElapsed)]);
%             [anglex, angley] = fx_decomSS_surface_roomsaving(ST_in, obj.front_I, obj.rear_I);
            [~, anglex, angley, anglez] = fx_decomSS_surface(ST_in, obj.front_I, obj.rear_I);
%             anglez = fx_decomSS_surface_zangle(ST_in, obj.front_I, obj.rear_I);
            % release the memory
            clear Gx;
            clear Gy;
            clear Gz;
            clear img;
            clear img_deconv;
            clear data;
            clear inam;
            % time record
            timeElapsed = toc;
            disp(['extraxt angles: ', num2str(timeElapsed)]);
            % interp3
%             obj.c_p     = cell2mat(cp);
            obj.angle_x = cell2mat(anglex);
            obj.angle_y = cell2mat(angley);
            obj.angle_z = cell2mat(anglez);
            % time record
            timeElapsed = toc;
            disp(['convert to arrays: ', num2str(timeElapsed)]);
        end
        
        function show_angles_ST(obj, medf_kernel, xslice, yslice, zslice, ds_rate)
            % display the extracted angles of ST
            % filter the noise by medfilter
            % medf_kernel: the kernel of median filter for smoothing
            % xslice, yslice, zslice: the slice in the 3d visualization;
            % ds_rate: the rate been applied for downsampling
            % these should be real unit here: mm, mm, us
            % for example:   xslice =  2; yslice = [2 5]; zslice = [];
            cp           = medfilt3(obj.c_p, medf_kernel, 'replicate');
            anglex       = medfilt3(real(obj.angle_x), medf_kernel, 'replicate');
            angley       = medfilt3(real(obj.angle_y), medf_kernel, 'replicate');
            [lx, ly, lz] = size(anglex);
            y            = (1: 1: lx*ds_rate) / obj.fx * 1e3;
            x            = (1: 1: ly*ds_rate) / obj.fy * 1e3;
            z            = (1: 1: lz) / obj.fs * 1e6;
            % ******************** show ply tracks
            cf    = figure;
            ax    = subplot(1, 1, 1);
            x_idx = xslice * obj.fx / 1e3;
            y_idx = yslice * obj.fy / 1e3;
            mask  = (obj.col==x_idx) | (obj.row==y_idx);
            x_row = obj.row(mask) / obj.fx * 1e3;
            x_col = obj.col(mask) / obj.fy * 1e3;  % / obj.fy * 1e3 ;
            x_dep = obj.dep(mask) / obj.fs * 1e6;  % / obj.fy * 1e3 ;
            scatter3(ax, x_col, x_row, x_dep, ...
                3,  [122 122 121]/255, 'filled', ...
                'DisplayName','Interply track');
            hold on;
            scatter3(xslice * ones(1,lx*ds_rate), ...
                y, obj.front_I(:, x_idx)/ obj.fs * 1e6, ...
                3, 'red', 'filled', ...
                'DisplayName','Front surface');
            hold on;
            scatter3(xslice * ones(1, lx*ds_rate), ...
                y, obj.rear_I(:, x_idx)/ obj.fs * 1e6, ...
                3, 'magenta', 'filled', ...
                'DisplayName','Rear surface');
            hold on;
            % select the yslice
            scatter3(x, yslice * ones(1, ly*ds_rate), ...
                obj.front_I(y_idx, :)/ obj.fs * 1e6, ...
                3, 'red', 'HandleVisibility','off');
            hold on;
            scatter3(x, yslice * ones(1, ly*ds_rate), ...
                obj.rear_I(y_idx, :)/ obj.fs * 1e6, ...
                3, 'magenta', 'HandleVisibility','off');
            hold on;
            % ****************************** anglex
            [X, Y, Z] = meshgrid(x(1:ds_rate:end), y(1:ds_rate:end), z);
            h         = slice(X, Y, Z, anglex, xslice, yslice, zslice);
            set(h,'EdgeColor','none');
            colormap jet;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle(degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            title('X-Z angle', 'fontsize', 16, 'Fontname', 'times new Roman');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            set(gca,'zDir','reverse'); % get time domain upside down
            lgd = legend('Interply', 'Front surface', 'Rear surface');
            % ************** show ply tracks
            cf    = figure;
            ax    = subplot(1, 1, 1);
            x_idx = xslice * obj.fx / 1e3;
            y_idx = yslice * obj.fy / 1e3;
            mask  = (obj.col==x_idx) | (obj.row==y_idx);
            x_row = obj.row(mask) / obj.fx * 1e3;
            x_col = obj.col(mask) / obj.fy * 1e3;  % / obj.fy * 1e3 ;
            x_dep = obj.dep(mask) / obj.fs * 1e6;  % / obj.fy * 1e3 ;
            scatter3(ax, x_col, x_row, x_dep, ...
                3,  [122 122 121]/255, 'filled', ...
                'DisplayName','Interply track');
            hold on;
            scatter3(xslice * ones(1, lx*ds_rate), ...
                y, obj.front_I(:, x_idx)/ obj.fs * 1e6, ...
                3, 'red', 'filled', ...
                'DisplayName','Front surface');
            hold on;
            scatter3(xslice * ones(1, lx*ds_rate), ...
                y, obj.rear_I(:, x_idx)/ obj.fs * 1e6, ...
                3, 'magenta', 'filled', ...
                'DisplayName','Rear surface');
            hold on;
            % select the yslice
            scatter3(x, yslice * ones(1, ly*ds_rate), ...
                obj.front_I(y_idx, :)/ obj.fs * 1e6, ...
                3, 'red', 'HandleVisibility','off');
            hold on;
            scatter3(x, yslice * ones(1, ly*ds_rate), ...
                obj.rear_I(y_idx, :)/ obj.fs * 1e6, ...
                3, 'magenta', 'HandleVisibility','off');
            hold on;
            % ************************ angley
            h = slice(X, Y, Z, angley, xslice, yslice, zslice);
            set(h,'EdgeColor','none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            title('Y-Z angle', 'fontsize', 16, 'Fontname', 'times new Roman');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            set(gca,'zDir','reverse'); % get time domain upside down
            lgd = legend('Interply', 'Front surface', 'Rear surface');
%             % c_p
%             cf = figure; h = slice(X, Y, Z, cp, xslice, yslice, zslice);
%             ax = subplot(1, 1, 1);
%             set(h,'EdgeColor','none');
%             colormap jet;
%             h = colorbar;
%             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle(degree)');
%             set(ax, 'fontsize', 16);
%             set(ax, 'linewidth', 1.5);
%             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
%             ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
%             zlabel('\fontname {times new roman} Time (\mus)', 'fontsize', 16);
%             title('C_p', 'fontsize', 16, 'Fontname', 'times new Roman');
%             set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
%             set(gca,'zDir','reverse'); % get time domain upside down
%             % caxis([-0.5, 0.5]);
        end
        
        function show_angles_ST_zangle(obj, medf_kernel, xslice, yslice, zslice, ds_rate)
            % display the extracted angles of ST
            % filter the noise by medfilter
            % medf_kernel: the kernel of median filter for smoothing
            % xslice, yslice, zslice: the slice in the 3d visualization;
            % ds_rate: the rate been applied for downsampling
            % these should be real unit here: mm, mm, us
            % for example:   xslice =  2; yslice = [2 5]; zslice = [];
            anglez       = medfilt3(real(obj.angle_z), medf_kernel, 'replicate');
            anglez       = (anglez + pi/2) / pi * 180;
            [lx, ly, lz] = size(anglez);
            y            = (1: 1: lx*ds_rate) / obj.fx * 1e3;
            x            = (1: 1: ly*ds_rate) / obj.fy * 1e3;
            z            = (1: 1: lz) / obj.fs * 1e3 * 3000/2;
            % ******************** show ply tracks
            cf = figure('Name', ['volumetric_zangles', '_', num2str(xslice)]);
            set(cf, 'Position', [0, 0, 1000, 800], 'color', 'white');
            ax    = subplot(1, 1, 1);
            x_idx = xslice * obj.fx / 1e3;
            y_idx = yslice * obj.fy / 1e3;
            mask  = (obj.col==x_idx) | (obj.row==y_idx);
%             x_row = obj.row(mask) / obj.fx * 1e3;
%             x_col = obj.col(mask) / obj.fy * 1e3;  % / obj.fy * 1e3 ;
%             x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
%             h1 = scatter3(ax, x_col, x_row, x_dep, ...
%                 3,  [122 122 121]/255, 'filled', ...
%                 'DisplayName','Interply track');
%             hold on;
            h2 = scatter3(xslice * ones(1,lx*ds_rate), ...
                y, obj.front_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
                3, 'red', 'filled', ...
                'DisplayName','Front surface');
            hold on;
            h3 = scatter3(xslice * ones(1, lx*ds_rate), ...
                y, obj.rear_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
                3, 'magenta', 'filled', ...
                'DisplayName','Rear surface');
            hold on;
            % select the yslice
            scatter3(x, yslice * ones(1, ly*ds_rate), ...
                obj.front_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
                3, 'red', 'HandleVisibility','off');
            hold on;
            scatter3(x, yslice * ones(1, ly*ds_rate), ...
                obj.rear_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
                3, 'magenta', 'HandleVisibility','off');
            hold on;
%             legend([h1 h2 h3], 'Interply track', 'Front surface', 'Rear surface');
            % ****************************** anglez
            [X, Y, Z] = meshgrid(x(1:ds_rate:end), y(1:ds_rate:end), z);
            h         = slice(X, Y, Z, anglez, xslice, yslice, zslice);
            set(h,'EdgeColor','none');
            colormap jet;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \alpha (\circ)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
            title('Out-of-plane angle', 'fontsize', 16, 'Fontname', 'times new Roman');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            set(gca,'zDir','reverse'); % get time domain upside down
            %             lgd = legend('Interply track', 'Front surface', 'Rear surface');
            legend([h2 h3], 'Front surface', 'Rear surface');
            xlim([0 x(end)]);
            ylim([0 y(end)]);
            zlim([1.5 8]);     
            zticks([1 3 5]);
            zticklabels({'0','2','4'});
            view([15 65 40]);
        end
        
        function show_angles_ST_Bscan(obj, B_type, index, win)
            % demonstrate the B_scan of the interply track
            % show inph_ex
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % win: the win of the B-scan image
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            x           = (1:(win(end)-win(1)+1))/obj.fx * 1e3;
            z           = (1:size(obj.img, 3))/obj.fs * 1e3 * 3000/2;
            cf          = figure('Name', ['B_scan_angles_', B_type, '_', num2str(index)]);
            set(cf, 'Position', [0, 0, length(win)*2, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            if (B_type == 'x')
                % B scan image    
                B_scan      = real(squeeze(obj.angle_x(index, win, :))) / pi *180;
%                 imagesc(x, z, B_scan');
                pcolor(x, z, B_scan');
                shading flat;
                set(gca, 'ydir', 'reverse');
                hold on;
                mask  = (obj.row==index) & obj.col>=win(1) & obj.col<=win(end);
                x_row = obj.col(mask) / obj.fx * 1e3;
                x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
%                 h1 = scatter(ax, x_row, x_dep, ...
%                     3,  [122 122 121]/255, 'filled', ...
%                     'DisplayName','Interply track');
%                 hold on;
                h2 = scatter(x, obj.front_I(index, win)/ obj.fs * 1e3 * 3000/2, ...
                    3, 'red', 'filled', ...
                    'DisplayName','Front surface');
                hold on;
                h3 = scatter(x, obj.rear_I(index, win)/ obj.fs * 1e3 * 3000/2, ...
                    3, 'magenta', 'filled', ...
                    'DisplayName','Rear surface');
                hold on;
                % figure setup
                colormap(jet);
                h = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (\circ)');
                %             cl   = caxis;
                %             caxis(cl);
                set(gca, 'Fontname', 'times new Roman');
                xlabel('\fontname {times new roman} x location (mm)', 'fontsize', 16);
                ylabel('\fontname {times new roman} z location (mm)', 'fontsize', 16);
                title('X-Z angle', 'fontsize', 16, 'Fontname', 'times new Roman');
            elseif (B_type == 'y')
                % B scan image
                B_scan      = real(squeeze(obj.angle_y(win, index, :))) / pi *180;
                %                 imagesc(x, z, B_scan');
                pcolor(x, z, B_scan');
                shading flat;
                set(gca, 'ydir', 'reverse');
                hold on;
                mask  = (obj.col==index) & obj.row>=win(1) & obj.row<=win(end);
                x_row = obj.row(mask) / obj.fy * 1e3;
                x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
%                 h1 = scatter(ax, x_row, x_dep, ...
%                     3,  [122 122 121]/255, 'filled', ...
%                     'DisplayName','Interply track');
%                 hold on;
                h2 = scatter(x, obj.front_I(win, index)/ obj.fs * 1e3 * 3000/2, ...
                    3, 'red', 'filled', ...
                    'DisplayName','Front surface');
                hold on;
                h3 = scatter(x, obj.rear_I(win, index)/ obj.fs * 1e3 * 3000/2, ...
                    3, 'magenta', 'filled', ...
                    'DisplayName','Rear surface');
                hold on;
                % figure setup
                colormap(jet);
                h = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (\circ)');
                %             cl   = caxis;
                %             caxis(cl);
                set(gca, 'Fontname', 'times new Roman');
                xlabel('\fontname {times new roman} y location (mm)', 'fontsize', 16);
                ylabel('\fontname {times new roman} z location (mm)', 'fontsize', 16);
                title('Y-Z angle', 'fontsize', 16, 'Fontname', 'times new Roman');
            end
            % legend([h1 h2 h3], 'Interply track', 'Front surface', 'Rear surface');
            legend([h2 h3], 'Interply track', 'Front surface', 'Rear surface');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 2);
            ylim([0.8 8]);
        end
        
        function makeMovie_angles_ST_Bscan(obj, B_type, win, property_name)
            % make a movie of the B_scan of the interply track
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % win: the win of the B-scan image
            % property_name: name of the dataset
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            dataset = obj.(property_name);
            dataset = medfilt3(dataset, [3, 3, 3], 'replicate');
            
            dataset = (dataset + pi/2) / pi * 180;
            front_p = obj.front_I;
            bakc_p  = obj.rear_I;
            x       = (1:(win(end)-win(1)+1))/obj.fx * 1e3;
            y       = (1:(win(end)-win(1)+1))/obj.fy * 1e3;
            z       = (1:size(dataset, 3))/obj.fs * 1e3 * 3000/2;
            cf      = figure('Name', ['B_scan_angles_', B_type, '_', num2str(win(end))]);
            set(cf, 'Position', [0, 0, length(win)*2, 1200], 'color', 'white');
            % make vedio
            axis tight manual;
            title(['y =' num2str(0/obj.fy * 1e3) ' mm']);
            xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} z (mm)', 'fontsize', 16);
            v = VideoWriter([property_name, '.avi']);
            open(v);
            %
            if (B_type == 'x')
                for index = 1:size(dataset,1)
                    % B scan image
                    % side view
                    B_scan = squeeze(dataset(index, win, :));
                    figure(cf);
                    ca1 = subplot(2, 1, 1);
                    % imagesc(x, z, B_scan');
                    pcolor(x, z, B_scan');
                    shading flat;
                    set(ca1, 'ydir', 'reverse');
                    hold on;
                    title(['y =' num2str(index/obj.fy * 1e3) ' mm']);
                    xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
                    ylabel('\fontname {times new roman} z (mm)', 'fontsize', 16);
                    %                     mask  = (obj.row==index) & obj.col>=win(1) & obj.col<=win(end);
                    %                     x_row = obj.col(mask) / obj.fx * 1e3;
                    %                     x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
                    %                     h1 = scatter(ca, x_row, x_dep, ...
                    %                         3,  [122 122 121]/255, 'filled', ...
                    %                         'DisplayName','Interply track');
                    %                     hold on;
                    h2 = scatter(x, front_p(index, win)/ obj.fs * 1e3 * 3000/2, ...
                        3, 'red', 'filled', ...
                        'DisplayName','Front surface');
                    hold on;
                    h3 = scatter(x, bakc_p(index, win)/ obj.fs * 1e3 * 3000/2, ...
                        3, 'magenta', 'filled', ...
                        'DisplayName','Rear surface');
                    hold on;
                    % figure setup
                    colormap(jet);
                    caxis(ca1, [0 4.5]);
                    h = colorbar;
                    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (\circ)');
                    set(ca1, 'Fontname', 'times new Roman');
                    set(ca1, 'Fontsize', 16);
                    set(ca1, 'Linewidth', 1.5);
                    % top view
                    figure(cf);
                    ca2     = subplot(2, 1, 2);
                    % imagesc(x, z, B_scan');
                    z_index     = round(mean(front_p, 'all', 'omitnan') ...
                        + mean(bakc_p-front_p, 'all', 'omitnan') / size(dataset,1) * index);
                    B_scan = squeeze(dataset(1:end, 1:end, z_index));
                    pcolor(x, y, B_scan);
                    shading flat;
                    hold on;
                    title(['z =' num2str(z_index/obj.fs * 1e3 * 3000/2) ' mm']);
                    xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
                    ylabel('\fontname {times new roman} y (mm)', 'fontsize', 16);
                    %                     mask  = (obj.row==index) & obj.col>=win(1) & obj.col<=win(end);
                    %                     x_row = obj.col(mask) / obj.fx * 1e3;
                    %                     x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
                    %                     h1 = scatter(ca, x_row, x_dep, ...
                    %                         3,  [122 122 121]/255, 'filled', ...
                    %                         'DisplayName','Interply track');
                    %                     hold on;
                    hold on;
                    % figure setup
                    colormap(jet);
                    caxis(ca2, [0 4.5]);
                    h = colorbar;
                    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (\circ)');
                    set(ca2, 'Fontname', 'times new Roman');
                    set(ca2, 'Fontsize', 16);
                    set(ca2, 'Linewidth', 1.5);
                    % frame vedio operation
                    frame = getframe(cf);
                    writeVideo(v, frame);
                    clf(cf);
                end
            elseif (B_type == 'y')
%                 % B scan image
%                 B_scan      = real(squeeze(obj.angle_y(win, index, :))) / pi *180;
%                 %                 imagesc(x, z, B_scan');
%                 pcolor(x, z, B_scan');
%                 shading flat;
%                 set(gca, 'ydir', 'reverse');
%                 hold on;
%                 mask  = (obj.col==index) & obj.row>=win(1) & obj.row<=win(end);
%                 x_row = obj.row(mask) / obj.fy * 1e3;
%                 x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
%                 h1 = scatter(ax, x_row, x_dep, ...
%                     3,  [122 122 121]/255, 'filled', ...
%                     'DisplayName','Interply track');
%                 hold on;
%                 h2 = scatter(x, obj.front_I(win, index)/ obj.fs * 1e3 * 3000/2, ...
%                     3, 'red', 'filled', ...
%                     'DisplayName','Front surface');
%                 hold on;
%                 h3 = scatter(x, obj.rear_I(win, index)/ obj.fs * 1e3 * 3000/2, ...
%                     3, 'magenta', 'filled', ...
%                     'DisplayName','Rear surface');
%                 hold on;
%                 % figure setup
%                 colormap(jet);
%                 h = colorbar;
%                 set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (\circ)');
%                 %             cl   = caxis;
%                 %             caxis(cl);
%                 set(gca, 'Fontname', 'times new Roman');
%                 xlabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
%                 ylabel('\fontname {times new roman} Z displacement (mm)', 'fontsize', 16);
%                 title('Y-Z angle', 'fontsize', 16, 'Fontname', 'times new Roman');
            end
%             legend([h1 h2 h3], 'Interply track', 'Front surface', 'Rear surface');
%             set(ca, 'fontsize', 16);
%             set(ca, 'linewidth', 2);
%             ylim([0.8 7]);
        end
        
        function obj = Filter_logGabor(obj, f0, sigma, PropertyName)
            % Filter the 3d dataset by logGabor filter
            % f0: the array of the center frequencies of the filter
            % sigma: affects the bandwidth of the filter.
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            img_ori = obj.(PropertyName);
            [lx, ly, n] = size(img_ori);
            % save in the obj.img_hil_filter
            img_temp = zeros(size(obj.(PropertyName)));
            % define a logGabor filter and calcualte the analytical signal
            fshift = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            filter =  fx_1dLogGabor_filter(fshift, f0, sigma);
            for x = 1: lx
                tic;
                for y = 1:ly
                    Ascan     = squeeze(img_ori(x, y, :));
                    % calculate the AS in frequency domain
                    Fo        = fft(Ascan);
                    % transverse Fo for dot multiply
                    Fo        = Fo';
                    fx_shift  = fftshift(Fo);
                    Fa        = (1 + sign(fshift)) .* fx_shift;
                    % add the logGabor filter here
                    Fa_filter = Fa .* filter;
                    % ifft to obtain the AS
                    Fa_ishift = ifftshift(Fa_filter);  % back to the original transform output
                    Ascan_as  = ifft(Fa_ishift);
                    % the direction of 'as' is reversed somehow, needs to be
                    % flipped and conjugated!
                    Ascan_as  = conj(flip(Ascan_as));
                    img_temp(x, y, :) = Ascan_as;   
                end
                clc;
                disp([num2str(x), '/' , num2str(lx)]);
            end
            obj.img_hil_filter = img_temp;
        end
        
        function obj = show_unwraped_inph_3d(obj, PropertyName, xslice, yslice, zslice, sigma1)
            % demonstrate the 3d  unwrapped phase
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % xslice =  100 / obj.fx * 1e3;
            % yslice = 50 / obj.fy * 1e3;
            % zslice = [];
            % sigma1: smoothing scale
            % ***
            % calculate the inph and unwrap it by each A-scan
            inph = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            inph_unw = zeros(size(inph));
            inph_unw_rear = zeros(size(inph(:, :, 1)));
            for idx = 1:size(inph, 1)
                for idy = 1:size(inph, 2)
                    z1 = round(obj.front_I(idx, idy));
                    %                     z2 = round(obj.rear_I(idx, idy)) + 50; % for smooth
                    A_inph = squeeze(inph(idx, idy, z1:end));
                    % track the phase by unwraped inst. phase
                    A_inph_unw = unwrap(A_inph);
                    inph_unw(idx, idy, z1: end) = A_inph_unw;
                end
            end
            y = (0: size(inph_unw, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_unw, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_unw, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            inph_visual = medfilt3(inph_unw, sigma1);
            % save it to the obj's property
            obj.inph_unwphase = inph_visual;
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure();
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            h = slice(ax, X, Y, Z, inph_visual, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. phase(rad.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            set(gca,'zDir','reverse'); % get time domain upside down
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            % the unwrapped phase of the rear surface
            for idx = 1:size(inph, 1)
                for idy = 1:size(inph, 2)
                    if ~isnan(obj.rear_I(idx, idy)) && obj.rear_I(idx, idy) < size(inph_visual, 3)
                        inph_unw_rear(idx, idy) = inph_visual(idx, idy, obj.rear_I(idx, idy));
                    end
                end
            end
            cf = figure();
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            imagesc(y, x, inph_unw_rear);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. phase(rad.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            set(gca,'zDir','reverse'); % get time domain upside down
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
        end
        
        function save_STangles_3D(obj, savename, supinfo)
            % save all the file of the 3D orientation
            % savename: the name incl. the path and the file for saving
            % save member variables
            front_I_var = obj.front_I;
            rear_I_var = obj.rear_I;
            angle_x_var = obj.angle_x;
            angle_y_var = obj.angle_y;
            img_hil_var = obj.img_hil;
            % supplement info.
            supinfo_var = supinfo;
            save(savename, ...
                'front_I_var', ...
                'rear_I_var', ...
                'angle_x_var', ...
                'angle_y_var', ...
                'img_hil_var', ...
                'supinfo_var', ...
                '-v7.3');
        end
        
        % *********************** ply track ***************
        function obj = track_interply(obj, PropertyName)
            % track the interply by unwrapped inph.
            % also derive the rear_I by phase
            % return a masked dataset
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % ************
            % The sequence is changed unexpectedly. Transverse the matrix
            inph               = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            img_temp           = obj.(PropertyName);
            mask_plytrack_temp = zeros(size(inph));
            row_temp           = [];
            col_temp           = [];
            dep_temp           = [];
            layer_counts_temp  = [];
            thickness          = NaN(size(mask_plytrack_temp,1), ...
                size(mask_plytrack_temp,2), 100); % unlimited
            % ************ filter the inph
%             3D Gaussian filter
%             sigma1             = [3 3 3];
%             img_temp_real      = real(img_temp);
%             img_temp_imag      = imag(img_temp);
%             img_temp_real_f    = imgaussfilt3(img_temp_real, sigma1, 'padding', 'replicate');
%             img_temp_imag_f    = imgaussfilt3(img_temp_imag, sigma1, 'padding', 'replicate');
%             img_temp           = img_temp_real_f + 1j*img_temp_imag_f;
            for idx = 1:size(inph, 1)
                for idy = 1:size(inph, 2)
                    if obj.rear_I(idx, idy) <=0 || obj.rear_I(idx, idy) > size(inph, 3) ...
                            || obj.front_I(idx, idy) <=0 || obj.front_I(idx, idy) > obj.rear_I(idx, idy)
                        continue;
                    end
                    front_I_interp                                   = obj.front_I(idx, idy);
                    rear_I_interp                                    = obj.rear_I(idx, idy);
                    Ascan_as                                         = squeeze(img_temp(idx, idy, :));
                    %             % track the - pi/2 by interplation on the unwraped data
                    Idx_resin                                        = linspace(front_I_interp, rear_I_interp, 26); % 25 for nominal value 
                    threshold                                        = 0.05;
                    [I_phase_track, ~, ~, ~, ~]                      = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);
                    obj.rear_I_pre(idx, idy)                         = obj.rear_I(idx, idy);
                    obj.rear_I(idx, idy)                             = I_phase_track(end);
                    I_phase_track                                    = I_phase_track(2:end-1);
                    mask_plytrack_temp(idx, idy, ...
                        round(I_phase_track(~isnan(I_phase_track)))) = 1;
                    % thickness calculation
                    thickness(idx, idy, 1:length(I_phase_track)-1)   = diff(I_phase_track);
                    % save the 3D scatter points
                    % updata the temp. context
                    layer_counts_temp = cat(2, layer_counts_temp, (1:length(I_phase_track)));
                    row_temp          = cat(2, row_temp,          idx*ones(1, length(I_phase_track)));
                    col_temp          = cat(2, col_temp,          idy*ones(1, length(I_phase_track)));
                    dep_temp          = cat(2, dep_temp,          round(I_phase_track));
                end
                clc;
                disp([num2str(idx), '/' , num2str(size(inph, 1))]);
            end
            obj.mask_plytrack = mask_plytrack_temp;
            obj.row           = row_temp;
            obj.col           = col_temp;
            obj.dep           = dep_temp;
            obj.layer_counts  = layer_counts_temp;
            obj.est_thickness = thickness;
        end
        
        function obj = track_interply_2ndharmonic(obj, PropertyName, nol)
            % track the interply by unwrapped inph.
            % also derive the rear_I by phase
            % return a masked dataset
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % nol: number of layers.
            % ************
            % The sequence is changed unexpectedly. Transverse the matrix
            inph               = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            img_temp           = obj.(PropertyName);
            mask_plytrack_temp = zeros(size(inph));
            row_temp           = [];
            col_temp           = [];
            dep_temp           = [];
            layer_counts_temp  = [];
            thickness          = NaN(size(mask_plytrack_temp,1), ...
                size(mask_plytrack_temp,2), 100); % unlimited
%             % ************ filter the inph
%             % 3D Gaussian filter
%             sigma1             = [5 5 5];
%             img_temp_real      = real(img_temp);
%             img_temp_imag      = imag(img_temp);
%             img_temp_real_f    = imgaussfilt3(img_temp_real, sigma1, 'padding', 'replicate');
%             img_temp_imag_f    = imgaussfilt3(img_temp_imag, sigma1, 'padding', 'replicate');
%             img_temp           = img_temp_real_f + 1j*img_temp_imag_f;
            %             inph_cos        = cos(inph);
%             inph_cos_filter = imgaussfilt3(inph_cos, sigma1, 'padding', 'replicate');
%             inph            = acos(inph_cos_filter);
            for idx = 1:size(inph, 1)
                for idy = 1:size(inph, 2)
                    if obj.rear_I(idx, idy) <=0 || obj.rear_I(idx, idy) > size(inph, 3) ...
                            || obj.front_I(idx, idy) <=0 || obj.front_I(idx, idy) > obj.rear_I(idx, idy)
                        continue;
                    end
                    % for debug
%                     if idx==22&&idy==100
%                         ;
%                     end
                    front_I_interp                                   = obj.front_I(idx, idy);
                    rear_I_interp                                    = obj.rear_I(idx, idy);
                    Ascan_as                                         = squeeze(img_temp(idx, idy, :));
                    %             % track the - pi/2 by interplation on the unwraped data
                    Idx_resin                                        = linspace(front_I_interp, rear_I_interp, nol);
                    threshold                                        = 0.05;
                    [I_phase_track, ~, ~, ~, ~]                      = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);
                    obj.rear_I_pre(idx, idy)                         = obj.rear_I(idx, idy);
                    obj.rear_I(idx, idy)                             = I_phase_track(end);
%                     I_phase_track                                    = I_phase_track(2:1:end-1);
                    I_phase_track                                    = I_phase_track(3:2:end-2); % 2nd harmonic, 2 tracks in one ply
                    mask_plytrack_temp(idx, idy, ...
                        round(I_phase_track(~isnan(I_phase_track)))) = 1;
                    % thickness calculation
                    thickness(idx, idy, 1:length(I_phase_track)-1)   = diff(I_phase_track);
                    % save the 3D scatter points
                    % updata the temp. context
                    layer_counts_temp = cat(2, layer_counts_temp, (1:length(I_phase_track)));
                    row_temp          = cat(2, row_temp,          idx*ones(1, length(I_phase_track)));
                    col_temp          = cat(2, col_temp,          idy*ones(1, length(I_phase_track)));
                    dep_temp          = cat(2, dep_temp,          round(I_phase_track));
                end
                disp([num2str(idx), '/' , num2str(size(inph, 1))]);
            end
            obj.mask_plytrack = mask_plytrack_temp;
            obj.row           = row_temp;
            obj.col           = col_temp;
            obj.dep           = dep_temp;
            obj.layer_counts  = layer_counts_temp;
            obj.est_thickness = thickness;
        end
        
        function obj = track_interply_hybrid(obj, PropertyName, f0_1, sigma0_1, f0_2, sigma0_2, nol)
            % track the interply by unwrapped inph.
            % also derive the rear_I by phase
            % return a masked dataset
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % nol: number of the interply layers
            % ************
            % The sequence is changed unexpectedly. Transverse the matrix
            img_temp           = obj.(PropertyName);
            [lx, ly, n]        = size(img_temp);
            mask_plytrack_temp = zeros([lx, ly, n]);
            row_temp           = [];
            col_temp           = [];
            dep_temp           = [];
            layer_counts_temp  = [];
            thickness          = NaN(size(mask_plytrack_temp,1), ...
                size(mask_plytrack_temp,2), 100); % unlimited
            % ************ filter the inph
            % 3D Gaussian filter
            sigma1             = [3 3 3];
            img_temp_real      = real(img_temp);
            img_temp_imag      = imag(img_temp);
            img_temp_real_f    = imgaussfilt3(img_temp_real, sigma1, 'padding', 'replicate');
            img_temp_imag_f    = imgaussfilt3(img_temp_imag, sigma1, 'padding', 'replicate');
            img_temp           = img_temp_real_f + 1j*img_temp_imag_f;
            % define the first logGabor filter
            fshift  = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            filter1 =  fx_1dLogGabor_filter(fshift, f0_1, sigma0_1);
            % define the second logGabor filter
            filter2 =  fx_1dLogGabor_filter(fshift, f0_2, sigma0_2);
            for idx = 1:lx
                for idy = 1:ly
                    if obj.rear_I(idx, idy) <=0 || obj.rear_I(idx, idy) > n || ...
                            obj.front_I(idx, idy) <=0 || obj.front_I(idx, idy) > obj.rear_I(idx, idy)
                        continue;
                    end
                    front_I_interp = round(obj.front_I(idx, idy));
                    rear_I_interp  = round(obj.rear_I(idx, idy));
                    Ascan          = squeeze(img_temp(idx, idy, :));
                    % ************* 2nd harmonic for the 1st and last interply ****
                    % calculate the AS in frequency domain
                    Fo         = fft(Ascan);
                    % transverse Fo for dot multiply
                    Fo         = Fo';
                    fx_shift   = fftshift(Fo);
                    Fa         = (1 + sign(fshift)) .* fx_shift;
                    % add the logGabor filter1 here
                    Fa_filter  = Fa .* filter1;
                    % ifft to obtain the AS
                    Fa_ishift  = ifftshift(Fa_filter);  % back to the original transform output
                    Ascan_as   = ifft(Fa_ishift);
                    % the direction of 'as' is reversed somehow, needs to be
                    % flipped and conjugated!
                    Ascan_as   = conj(flip(Ascan_as));
                    % **** track the - pi/2 by interplation on the unwraped data ****
                    inph_ascan = angle(Ascan_as);
                    tol        = pi;
                    inph_ascan_unwrap = unwrap(inph_ascan, tol);
                    if rear_I_interp > 0 && front_I_interp > 0 && ...
                            inph_ascan_unwrap(rear_I_interp)-inph_ascan_unwrap(front_I_interp) >= 3 * pi
                        phase_seq     = inph_ascan_unwrap(round(front_I_interp))+pi: 2*pi: inph_ascan_unwrap(round(rear_I_interp))-pi; % back_inph - pi!
                        % remove duplicates using 'unique' function.
                        [~, ind]      = unique(inph_ascan_unwrap); % ind = index of first occurrence of a repeated value
                        v             = 1:length(inph_ascan_unwrap);
                        I_phase_track = interp1(inph_ascan_unwrap(ind), v(ind), phase_seq);
                        % ******** define the phase_0 - pi/2 ********
                        phase_0       = interp1(round(front_I_interp)+(-1:1), ...
                            inph_ascan(round(front_I_interp)+(-1:1)), front_I_interp);
                        % phase_0                 = inph_ascan(front_I);
                        % define phase_0 - pi/2
                        if phase_0-pi/2 >= -pi
                            phase_pi2 = phase_0 - pi/2;
                        else
                            phase_pi2 = phase_0 + pi*3/2;
                        end
                        lct = NaN(1, length(I_phase_track));
                        for n_layer = [2 length(I_phase_track)-1]
                            % reaches the true number of the layers, which means this layer should
                            if n_layer+1 > length(I_phase_track) || n_layer+1 < 2
                                continue;
                            end
                            idx_layer_n1 = round(I_phase_track(n_layer));
                            idx_layer_n2 = round(I_phase_track(n_layer+1));
                            % find one nearest point and apply interpolation
                            inph_ascan_oneply         = inph_ascan(idx_layer_n1: idx_layer_n2);
                            %in order to only find positive numbers
                            inph_check                = inph_ascan_oneply - phase_pi2;
                            inph_check(inph_check<=0) = nan; %replace negative numbers and the zero with nan
                            [~, I_minpos]             = min(inph_check); %find values
                            %in order to only find negative numbers
                            inph_check                = inph_ascan_oneply - phase_pi2;
                            inph_check(inph_check>0)  = nan; %replace negative numbers and the zero with nan
                            [~, I_maxneg]             = max(inph_check); %find values
                            nearest_ps                = [I_maxneg I_minpos];
                            x                         = inph_ascan_oneply(nearest_ps);
                            [~, ind]                  = unique(x); % ind = index of first occurrence of a repeated value
                            if isempty(ind)
                                lct(n_layer) = NaN;
                            elseif length(ind)==1 || ...
                                    isempty(find(x>phase_pi2, 1)) || ...
                                    isempty(find(x<=phase_pi2, 1))
                                lct(n_layer) = round(mean(nearest_ps(ind)))+idx_layer_n1-1;
                            else
                                vq1                  = interp1(x(ind), nearest_ps(ind), phase_pi2, 'linear');
                                lct(n_layer)         = vq1+idx_layer_n1-1;
                            end
                        end
                    else
                        continue;
                    end
                    if length(lct) > 2
                        interply_1st      = lct(2);
                        interply_last     = lct(end-1);
                    end
                    % ************* fundamental harmonic for the rest interplies ****
                    % calculate the AS in frequency domain
                    Fo       = fft(Ascan);
                    % transverse Fo for dot multiply
                    Fo        = Fo';
                    fx_shift  = fftshift(Fo);
                    Fa        = (1 + sign(fshift)) .* fx_shift;
                    % add the logGabor filter2 here
                    Fa_filter = Fa .* filter2;
                    % ifft to obtain the AS
                    Fa_ishift = ifftshift(Fa_filter);  % back to the original transform output
                    Ascan_as  = ifft(Fa_ishift);
                    % the direction of 'as' is reversed somehow, needs to be
                    % flipped and conjugated!
                    Ascan_as  = conj(flip(Ascan_as));
                    %                     lct       = NaN(1, nol);
                    if interply_last>interply_1st
                        Idx_resin = linspace(front_I_interp, rear_I_interp, nol);
                        threshold = 0.05;
                        [I_phase_track, ~, ~, ~, ~]         = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);
                        I_phase_track(isnan(I_phase_track)) = [];
                        I_phase_track = cat(2, interply_1st, I_phase_track(3:end-2), interply_last);
                    else % interply_last<=interply_1st only one interply
                        I_phase_track = interply_1st;
                    end
                    mask_plytrack_temp(idx, idy, ...
                        round(I_phase_track(~isnan(I_phase_track)))) = 1;
                    % thickness calculation
                    thickness(idx, idy, 1:length(I_phase_track)-1)   = diff(I_phase_track);
                    % save the 3D scatter points
                    % updata the temp. context
                    layer_counts_temp = cat(2, layer_counts_temp, (1:length(I_phase_track)));
                    row_temp          = cat(2, row_temp,          idx*ones(1, length(I_phase_track)));
                    col_temp          = cat(2, col_temp,          idy*ones(1, length(I_phase_track)));
                    dep_temp          = cat(2, dep_temp,          round(I_phase_track));
                end
                clc;
                disp([num2str(idx), '/' , num2str(lx)]);
            end
            obj.mask_plytrack = mask_plytrack_temp;
            obj.row           = row_temp;
            obj.col           = col_temp;
            obj.dep           = dep_temp;
            obj.layer_counts  = layer_counts_temp;
            obj.est_thickness = thickness;
        end
        
        function obj = track_interply_hybrid_ascan(obj, PropertyName, f0_1, sigma0_1, f0_2, sigma0_2, nol, idx, idy)
            % track the interply by unwrapped inph.
            % also derive the rear_I by phase
            % return a masked dataset
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % nol: number of the interply layers
            % idx, idy: index to select ascan
            % ************
            % The sequence is changed unexpectedly. Transverse the matrix
            img_temp         = obj.(PropertyName);
            [~, ~, n]        = size(img_temp);
            t            = (1:n);
            t_domain     = t / obj.fs*1e6;
            % ************ filter the inph
            % 3D Gaussian filter
            sigma1             = [3 3 3];
            img_temp_real      = real(img_temp);
            img_temp_imag      = imag(img_temp);
            img_temp_real_f    = imgaussfilt3(img_temp_real, sigma1, 'padding', 'replicate');
            img_temp_imag_f    = imgaussfilt3(img_temp_imag, sigma1, 'padding', 'replicate');
            img_temp           = img_temp_real_f + 1j*img_temp_imag_f;
            % define the first logGabor filter
            fshift  = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            filter1 =  fx_1dLogGabor_filter(fshift, f0_1, sigma0_1);
            % define the second logGabor filter
            filter2 =  fx_1dLogGabor_filter(fshift, f0_2, sigma0_2);
            front_I_interp = round(obj.front_I(idx, idy));
            rear_I_interp  = round(obj.rear_I(idx, idy));
            Ascan          = squeeze(img_temp(idx, idy, :));
            % ************* 2nd harmonic for the 1st and last interply ****
            % calculate the AS in frequency domain
            Fo         = fft(Ascan);
            % transverse Fo for dot multiply
            Fo         = Fo';
            fx_shift   = fftshift(Fo);
            Fa         = (1 + sign(fshift)) .* fx_shift;
            % add the logGabor filter1 here
            Fa_filter  = Fa .* filter1;
            % ifft to obtain the AS
            Fa_ishift  = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as   = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as   = conj(flip(Ascan_as));
            % **** track the - pi/2 by interplation on the unwraped data ****
            inph_ascan = angle(Ascan_as);
            tol        = pi;
            inph_ascan_unwrap = unwrap(inph_ascan, tol);
            if rear_I_interp > 0 && front_I_interp > 0 && ...
                    inph_ascan_unwrap(rear_I_interp)-inph_ascan_unwrap(front_I_interp) >= 3 * pi
                phase_seq     = inph_ascan_unwrap(round(front_I_interp))+pi: 2*pi: inph_ascan_unwrap(round(rear_I_interp))-pi; % back_inph - pi!
                % remove duplicates using 'unique' function.
                [~, ind]      = unique(inph_ascan_unwrap); % ind = index of first occurrence of a repeated value
                v             = 1:length(inph_ascan_unwrap);
                I_phase_track = interp1(inph_ascan_unwrap(ind), v(ind), phase_seq);
                % ******** define the phase_0 - pi/2 ********
                phase_0       = interp1(round(front_I_interp)+(-1:1), ...
                    inph_ascan(round(front_I_interp)+(-1:1)), front_I_interp);
                % phase_0                 = inph_ascan(front_I);
                % define phase_0 - pi/2
                if phase_0-pi/2 >= -pi
                    phase_pi2 = phase_0 - pi/2;
                else
                    phase_pi2 = phase_0 + pi*3/2;
                end
                lct = NaN(1, length(I_phase_track));
                for n_layer = [2 length(I_phase_track)-1]
                    % reaches the true number of the layers, which means this layer should
                    if n_layer+1 > length(I_phase_track) || n_layer+1 < 2
                        continue;
                    end
                    idx_layer_n1 = round(I_phase_track(n_layer));
                    idx_layer_n2 = round(I_phase_track(n_layer+1));
                    % find one nearest point and apply interpolation
                    inph_ascan_oneply         = inph_ascan(idx_layer_n1: idx_layer_n2);
                    %in order to only find positive numbers
                    inph_check                = inph_ascan_oneply - phase_pi2;
                    inph_check(inph_check<=0) = nan; %replace negative numbers and the zero with nan
                    [~, I_minpos]             = min(inph_check); %find values
                    %in order to only find negative numbers
                    inph_check                = inph_ascan_oneply - phase_pi2;
                    inph_check(inph_check>0)  = nan; %replace negative numbers and the zero with nan
                    [~, I_maxneg]             = max(inph_check); %find values
                    nearest_ps                = [I_maxneg I_minpos];
                    x                         = inph_ascan_oneply(nearest_ps);
                    [~, ind]                  = unique(x); % ind = index of first occurrence of a repeated value
                    if isempty(ind)
                        lct(n_layer) = NaN;
                    elseif length(ind)==1 || ...
                            isempty(find(x>phase_pi2, 1)) || ...
                            isempty(find(x<=phase_pi2, 1))
                        lct(n_layer) = round(mean(nearest_ps(ind)))+idx_layer_n1-1;
                    else
                        vq1                  = interp1(x(ind), nearest_ps(ind), phase_pi2, 'linear');
                        lct(n_layer)         = vq1+idx_layer_n1-1;
                    end
                end
            end
            if length(lct) > 2
                interply_1st      = lct(2);
                interply_last     = lct(end-1);
            end
            % plot in fig.
            cf = figure('Name', ['Ascan_', 'hybridphase_xy_', num2str(idx), '_', num2str(idy)]);
            set(cf, 'Position', [0, 0, 800, 900], 'color', 'white');
            ca = subplot(3, 1, 1);
            h1 = plot(ca, t_domain, Ascan, 'k-', 'linewidth', 2);
            hold on;
            h2 = plot(ca, t_domain, abs(Ascan), '-', 'linewidth', 2, 'Color', [1 0.5 0]); % orange color
            hold on;
            h3 = plot(ca, t_domain(front_I_interp), abs(Ascan(front_I_interp)), 'ro', 'linewidth', 2);
            hold on;
            plot(t_domain(rear_I_interp), abs(Ascan(rear_I_interp)), 'ro', 'linewidth', 2);
            hold on;        
            xlabel('\fontname {times new roman} Time (\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Amp. (arb.)', 'fontsize', 16);
            set(ca, 'fontsize', 16);  
            set(ca, 'fontname', 'Times new roman');
            set(ca, 'linewidth', 1.5);
            legend([h1 h2 h3], 'Original signal', 'Instantaneous amplitude', 'Surfaces');
            xlim([t_domain(1) t_domain(end)]);
            ca = subplot(3, 1, 2);
            h1 = plot(ca, t_domain, inph_ascan, 'b-', 'linewidth', 2);
            hold on;
            h2 = plot(ca, t_domain(front_I_interp), inph_ascan(front_I_interp), 'ro', 'linewidth', 2);
            hold on;
            plot(t_domain(rear_I_interp), inph_ascan(rear_I_interp), 'ro', 'linewidth', 2);
            hold on;
            h3 = plot(ca, t_domain(round([interply_1st interply_last])), ...
                inph_ascan(round([interply_1st interply_last])), 'yo', 'linewidth', 2);
            hold on;
            xlabel('\fontname {times new roman} Time (\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Amp. (arb.)', 'fontsize', 16);
            set(ca, 'fontsize', 16);  
            set(ca, 'fontname', 'Times new roman');
            set(ca, 'linewidth', 1.5);
            legend([h1 h2 h3], 'Instantaneous phase',  'Surfaces', '1st and last interplies');
            xlabel('\fontname {times new roman} Time (\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Phase (rad.)', 'fontsize', 16);
            set(ca, 'fontsize', 16);  
            set(ca, 'fontname', 'Times new roman');
            set(ca, 'linewidth', 1.5);
            xlim([t_domain(1) t_domain(end)]);
            % ************* fundamental harmonic for the rest interplies ****
            % calculate the AS in frequency domain
            Fo       = fft(Ascan);
            % transverse Fo for dot multiply
            Fo        = Fo';
            fx_shift  = fftshift(Fo);
            Fa        = (1 + sign(fshift)) .* fx_shift;
            % add the logGabor filter2 here
            Fa_filter = Fa .* filter2;
            % ifft to obtain the AS
            Fa_ishift = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as  = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as  = conj(flip(Ascan_as));
            %                     lct       = NaN(1, nol);
            if interply_last>interply_1st
                Idx_resin = linspace(front_I_interp, rear_I_interp, nol);
                threshold = 0.05;
                [I_phase_track, ~, ~, ~, ~]         = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);
                I_phase_track(isnan(I_phase_track)) = [];
                I_phase_track = cat(2, interply_1st, I_phase_track(3:end-2), interply_last);
            else % interply_last<=interply_1st only one interply
                I_phase_track = interply_1st;
            end
             ca = subplot(3, 1, 3);
            h1 = plot(ca, t_domain, angle(Ascan_as), 'b-', 'linewidth', 2);
            hold on;
            h2 = plot(ca, t_domain(front_I_interp), inph_ascan(front_I_interp), 'ro', 'linewidth', 2);
            hold on;
            plot(t_domain(rear_I_interp), inph_ascan(rear_I_interp), 'ro', 'linewidth', 2);
            hold on;
            h3 = plot(ca, t_domain(round([interply_1st interply_last])), ...
                inph_ascan(round([interply_1st interply_last])), 'yo', 'linewidth', 2);
            hold on;
            h4 = plot(ca, t_domain(round(I_phase_track(2:end-1))), ...
                angle(Ascan_as(round(I_phase_track(2:end-1)))), 'go', 'linewidth', 2);
            xlabel('\fontname {times new roman} Time (\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Amp. (arb.)', 'fontsize', 16);
            set(ca, 'fontsize', 16);  
            set(ca, 'fontname', 'Times new roman');
            set(ca, 'linewidth', 1.5);
            legend([h1 h2 h3], 'Instantaneous phase',  'Surfaces', '1st and last interplies');
            xlabel('\fontname {times new roman} Time (\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Phase (rad.)', 'fontsize', 16);
            set(ca, 'fontsize', 16);  
            set(ca, 'fontname', 'Times new roman');
            set(ca, 'linewidth', 1.5);
            xlim([t_domain(1) t_domain(end)]);
        end
        
        function obj = track_interply_inph(obj, threshold, PropertyName, nol)
            % track the interply by instantaneous phases.
            % increase step by one-ply TOF
            % also derive the rear_I by phase
            % return a masked dataset
            % threshold: the threshold to select the peaks from inam
            % phase close to -pi/2
            % PropertyName: 'img_WienerDeconv_B', 'img_WienerDeconv_AR_B',
            % or 'img_SparseDeconv_B ...
            % nol: number of layers
            % ************
            % The sequence is changed unexpectedly. Transverse the matrix
            img_ori            = obj.(PropertyName);
            mask_plytrack_temp = zeros(size(img_ori));
            lct_3d             = NaN(size(img_ori, 1), size(img_ori, 2), nol);                
            row_temp           = [];
            col_temp           = [];
            dep_temp           = [];
            layer_counts_temp  = [];
            errors_pos_repeat  = NaN(size(mask_plytrack_temp,1), ...
                size(mask_plytrack_temp,2), nol);
            thickness          = NaN(size(mask_plytrack_temp,1), ...
                size(mask_plytrack_temp,2), nol-1);
            for idx = 1:size(img_ori, 1)
                for idy = 1:size(img_ori, 2)
                    if obj.rear_I(idx, idy) <=0 || obj.rear_I(idx, idy) > size(img_ori, 3) ...
                            || obj.front_I(idx, idy) <=0 || obj.front_I(idx, idy) > obj.rear_I(idx, idy)
                        continue;
                    end
%                     if idx ==9 && idy==80 % for debugging
%                         a = 1;
%                     end
                    A_ori                  = squeeze(img_ori(idx, idy, :));
                    % time-varying gain
%                     A_ori                  = A_ori .* (1:length(A_ori))';
                    % track the phase by wrapped inst. phase
                    %                     Idx_resin                                  = linspace(obj.front_I(idx, idy), obj.rear_I(idx, idy), 25);
                    %                     [lct, ~, ~, ~, ~]                          = fx_track_byphase(Idx_resin, A_ori, threshold);
                    % track the -pi/2 by distance
                    %                     [lct, ~, ~] = fx_track_bycluster(hilbert(A_ori), threshold, 23);
                    Idx_resin              = linspace(obj.front_I(idx, idy), ...
                       obj.rear_I(idx, idy), 25);  % !!!!! 25 is the number of layers incl. front and back surfaces
                    % track the -pi/2 by phase
                    % [lct, ~, ~, ~, ~]      = fx_track_byphase(Idx_resin, A_ori, threshold);
                    power_noise            = 1;
                    [lct, pct, err, ~, ~]  = fx_track_byphase_step(Idx_resin, real(A_ori), obj.fs, threshold, power_noise);
                    errors_pos_repeat(idx, idy, :) = err;
                    % thickness calculation
                    thickness(idx, idy, :) = diff(lct);
                    % save to 3d vector
                    if sum(isnan(lct)) > 10
                        temp = 1;
                    end
                    lct_3d(idx, idy, :)    = round(lct);
%                     mask_plytrack_temp(idx, idy, lct(2:end-1)) = 1;
                    % derive the rear surface by phase
                    %                     I_rear_surface       = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                    %                         phase_rear_surface) + round(obj.front_I(idx, idy));
                    obj.rear_I(idx, idy)   = lct(end);
                    obj.front_I(idx, idy)  = lct(1);
                    % save the 3D scatter points
                    num_interplies         = length(lct)-2; % incl. the front and rear surfaces
                    layer_counts_temp      = cat(2, layer_counts_temp, (1:num_interplies));
                    row_temp               = cat(2, row_temp, idx * ones(1, num_interplies));
                    col_temp               = cat(2, col_temp, idy * ones(1, num_interplies));
                    dep_temp               = cat(2, dep_temp, round(lct(2:end-1)));
                end
                disp([num2str(idx), '/' , num2str(size(img_ori, 1))]);
            end
            lct_3d_interp = lct_3d;
%             % interpolation
%             for i = 1:25
%                 lct_3d_interp(:, :, i) = round(fx_inpaint_nans(lct_3d(:, :, i), 5));
%             end
            % get the ply thickness
            P             = polyfit(1:length(Idx_resin), Idx_resin, 1);
            Idx_oneply    = P(1);
            lct_3d_thick  = diff(lct_3d, 1, 3) - Idx_oneply;
            error_mean    = squeeze(mean(lct_3d_thick,[1 2], 'omitnan'));
            error_std     = squeeze(std(lct_3d_thick, 0, [1 2], 'omitnan'));
            figure, errorbar(error_mean*100, error_std*100,'-s',...
                'MarkerSize', 10, 'Linewidth', 3,'MarkerEdgeColor','red','MarkerFaceColor','red');  
            % by individual saved err
            error_mean    = squeeze(mean(abs(errors_pos_repeat),[1 2], 'omitnan'));
            error_std     = squeeze(std(errors_pos_repeat, 0, [1 2], 'omitnan'));
            cf            = figure('Name', ['error_bar' '_' PropertyName]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca            = subplot(1, 1, 1);
            errorbar((1:nol)/obj.fs*1e6, error_mean*100, error_std*100,'-s',...
                'MarkerSize', 10, 'Linewidth', 3,'MarkerEdgeColor','red','MarkerFaceColor','red')
            xlabel('\fontname{times new roman} Time (µm)', 'fontsize', 16);
            ylabel('\fontname{times new roman} Errors (mean and std)', 'fontsize', 16);
            ytickformat('percentage');
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            grid on;
            % fill mask vector
            for idx = 1:size(img_ori, 1)
                for idy = 1:size(img_ori, 2)
                    % ignore NaN
                    for idz = 2:size(lct_3d_interp, 3) - 1
                        if ~isnan(lct_3d_interp(idx, idy, idz))
                            mask_plytrack_temp(idx, idy, lct_3d_interp(idx, idy, idz)) = 1;
                        end
                    end
                end
            end
            obj.mask_plytrack = mask_plytrack_temp;
            obj.errors_pos    = errors_pos_repeat;
            obj.row           = row_temp;
            obj.col           = col_temp;
            obj.dep           = dep_temp;
            obj.layer_counts  = layer_counts_temp;
            obj.est_thickness = thickness;
        end
        
        function demo_Ascan_plytrack(obj, x, y, interp_factor)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % interp_factor: the interpolation multiple
            Ascan          = squeeze(obj.img(x, y, :));
            % interpolation
            xp             = 1: 1: length(Ascan);
            xq             = 1: 1 / interp_factor: length(Ascan);
            Ascan_interp   = interp1(xp, Ascan, xq);
            front_I_interp = obj.front_I(x, y) * interp_factor;
            rear_I_interp  = obj.rear_I(x, y) * interp_factor;
            % plot the inam and the record signal
            Ascan_as       = hilbert(Ascan_interp);
            cf             = figure('Name', strcat('plytrack', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca             = subplot(2, 1, 1);
            plot(ca, (1:length(Ascan_interp))/obj.fs*1e6, Ascan_interp, 'linewidth', 2);
            hold on;
            plot(ca, (1:length(Ascan_as))/obj.fs*1e6, abs(Ascan_as), 'linewidth', 2)
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            legend({'Recorded', 'Inst. Amp.'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % % track the - pi / 2 by interplation on the unwraped data
            % A_inph        = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            % A_inph_unw    = unwrap(A_inph);
            % phase_seq     = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
            % I_phase_track = interp1(A_inph_unw, 1:length(A_inph_unw), ...
            %   phase_seq) + front_I_interp;
            % track the - pi / 2 by interplation on the wraped data
            Idx_resin           = linspace(front_I_interp, rear_I_interp, 25);
            threshold           = 0.05;
            [lct, pct, ~, ~, ~] = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);     
            % drop the Nan in 'I_phase_track'
            I_phase_track       = lct(2:end-1);
            I_phase_track       = I_phase_track(~isnan(I_phase_track));
            P_phase_track       = pct(2:end-1);
            P_phase_track       = P_phase_track(~isnan(P_phase_track));
            % plot the inph, and the ply track of - pi / 2
            ca                  = subplot(2, 1, 2);
            plot(ca, (1:length(Ascan_as))/obj.fs*1e6, angle(Ascan_as), 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            % used extracted location of the back-surface here
            plot(ca, [front_I_interp lct(end)]/obj.fs*1e6 , angle(Ascan_as(round([front_I_interp lct(end)]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            hold on;
            plot(ca, I_phase_track/obj.fs*1e6, P_phase_track, 'gd', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16} Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            legend({'Inst. phase', 'Derived surface', 'Derived inter-ply'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function demo_Ascan_plytrack_inph(obj, x, y, interp_factor, threshold)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % interp_factor: the interpolation multiple
            % threshold: threshold is a ratio to find the peaks
            % 
            Ascan                  = squeeze(obj.img(x, y, :));
            % interpolation
            xp                     = 1: 1: length(Ascan);
            xq                     = 1: 1 / interp_factor: length(Ascan);
            Ascan_interp           = interp1(xp, Ascan, xq);
            % find the front and back
            inam_ascan             = abs(hilbert(Ascan_interp));
            MinPeakHeight          = threshold * max(inam_ascan);
            MinPeakDistance        = length(inam_ascan) / 3;
            [~, Idx_walls]         = findpeaks(inam_ascan, 1:length(inam_ascan), ...
                'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
            front_I_interp         = Idx_walls(1);
            rear_I_interp          = Idx_walls(end);
            % plot the inam and the record signal
            Ascan_as               = hilbert(Ascan_interp);
            %             % track the - pi/2 by interplation on the unwraped data
            Idx_resin              = linspace(front_I_interp, rear_I_interp, 25);
            %             [lct, ~, ~, ~, A_inph] = fx_track_byphase(Idx_resin, real(Ascan_as), threshold);
            % track the -pi/2 by distance
            %             [lct, ~, ~]          = fx_track_bycluster(A_ori, threshold, 23, Idx_resin(2:end-1));
            %             [lct, inam_ascan, Ascan_phase, errors_pos] = fx_track_bydistance(Ascan_as, threshold, Idx_resin);
            % track the -pi/2 by phase
            [lct, pct, ~, inam_ascan, Ascan_phase] = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);
            % plot the inph, and the ply track of - pi / 2
            cf             = figure('Name', strcat('plytrack', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca             = subplot(2, 1, 1);
            t_space        = (1:length(Ascan_as))/obj.fs*1e6;
            plot(ca, t_space, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, t_space, abs(Ascan_as), 'linewidth', 2)
            hold on;
            plot(ca, t_space([front_I_interp rear_I_interp]), ...
                inam_ascan([front_I_interp rear_I_interp]), 'ro', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            legend({'Recorded signal', 'Inst. Amp.', 'Determined positions of surfaces'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            ca = subplot(2, 1, 2);
            plot(ca, t_space, Ascan_phase, 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, lct([1 end])/obj.fs*1e6, pct([1 end]), 'ro', 'linewidth', 1);
            hold on;
            plot(ca, lct(2:end-1)/obj.fs*1e6, pct(2:end-1), 'ms', 'linewidth', 1);
            hold on;
            plot(ca, t_space(round(Idx_resin(2:end-1))), Ascan_phase(round(Idx_resin(2:end-1))), 'gd', 'linewidth', 1);
            hold on;
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            legend({'Inst. phase', 'Determined positions of surfaces', 'Determined positions of interplies', 'Theoretical positions of interplies'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
%             % for debug
%             figure,
%             plot(errors_pos);
        end
        
        function demo_Ascan_plytrack_inph_v2(obj, x, y, interp_factor, threshold)
            % demonstrate the A scan and the ply track with a line showing
            % -pi/w
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % interp_factor: the interpolation multiple
            % threshold: threshold is a ratio to find the peaks
            % 
            Ascan                  = squeeze(obj.img(x, y, :));
            % interpolation
            xp                     = 1: 1: length(Ascan);
            xq                     = 1: 1 / interp_factor: length(Ascan);
            Ascan_interp           = interp1(xp, Ascan, xq);
            % find the front and back
            inam_ascan             = abs(hilbert(Ascan_interp));
            MinPeakHeight          = threshold * max(inam_ascan);
            MinPeakDistance        = length(inam_ascan) / 2;
            [~, Idx_walls]         = findpeaks(inam_ascan, 1:length(inam_ascan), ...
                'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
            front_I_interp         = round(Idx_walls(1));
            rear_I_interp          = round(Idx_walls(end));
            % plot the inam and the record signal
            Ascan_as               = hilbert(Ascan_interp);
            %             % track the - pi/2 by interplation on the unwraped data
            Idx_resin              = linspace(front_I_interp, rear_I_interp, 25);
            %             [lct, ~, ~, ~, A_inph] = fx_track_byphase(Idx_resin, real(Ascan_as), threshold);
            % track the -pi/2 by distance
            %             [lct, ~, ~]          = fx_track_bycluster(A_ori, threshold, 23, Idx_resin(2:end-1));
            %             [lct, inam_ascan, Ascan_phase, errors_pos] = fx_track_bydistance(Ascan_as, threshold, Idx_resin);
            % track the -pi/2 by phase
            [lct, ~, errors_pos, inam_ascan, Ascan_phase] = fx_track_byphase(Idx_resin, Ascan_as, threshold);
            % plot the inph, and the ply track of - pi / 2
            cf             = figure('Name', strcat('plytrack_v2_', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca             = subplot(2, 1, 1);
            t_space        = (1:length(Ascan_as)) / obj.fs * 1e6;
            plot(ca, t_space, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, t_space, abs(Ascan_as), 'linewidth', 2)
            hold on;
            plot(ca, t_space([front_I_interp rear_I_interp]), ...
                inam_ascan([front_I_interp rear_I_interp]), 'ro', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Recorded signal', 'Inst. Amp.', 'Determined positions of surfaces'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            ca = subplot(2, 1, 2);
            plot(ca, t_space, Ascan_phase, 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, t_space(lct([1 end])), Ascan_phase(lct([1 end])), 'ro', 'linewidth', 1);
            hold on;
            plot(ca, t_space(round(Idx_resin(2:end-1))), Ascan_phase(round(Idx_resin(2:end-1))), 'gd', 'linewidth', 1);
            hold on;     
            yl = yline(ca, Ascan_phase(lct(1))-pi/2,'--','LineWidth',2);
            yl.Color = [.80 0 .40];
            text(t_space(end), Ascan_phase(lct(1))-pi/2+pi/5, '\phi_0-\pi/2', 'Color', [.80 0 .40],'FontSize',16);
            hold on;
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Inst. phase', 'Determined positions of surfaces', 'Theoretical positions of interplies'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % for debug
            figure,
            plot(errors_pos);
        end
        
        function demo_Ascan_plytrack_addnoise(obj, x, y, interp_factor, snr_value)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % interp_factor: the interpolation multiple
            % snr: snr of the added noise
            Ascan = squeeze(obj.img(x, y, :));
            % interpolation
            xp             = 1: 1: length(Ascan);
            xq             = 1: 1/interp_factor: length(Ascan);
            Ascan_interp   = interp1(xp, Ascan, xq);
            front_I_interp = obj.front_I(x, y) * interp_factor;
            rear_I_interp  = obj.rear_I(x, y) * interp_factor;
            % add noise
            Ascan_interp   = awgn(Ascan_interp, snr_value, 'measured');
            % calcualte the snr
            %             SNR = snr(real(Ascan_interp), obj.fs);
            pn       = bandpower(Ascan_interp, obj.fs, [10e6, round(obj.fs/2.1)]);
            ps       = bandpower(Ascan_interp, obj.fs, [1e6, 10e6]);
            SNR      = 10*log10(ps/pn);
            disp("SNR:");
            disp(SNR);
            % % track the - pi / 2 by interplation on the unwraped data
            % A_inph        = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            % A_inph_unw    = unwrap(A_inph);
            % phase_seq     = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
            % I_phase_track = interp1(A_inph_unw, 1:length(A_inph_unw), ...
            %   phase_seq) + front_I_interp;
            % track the - pi / 2 by interplation on the wraped data
            Ascan_as            = hilbert(Ascan_interp);
            Idx_resin           = linspace(front_I_interp, rear_I_interp, 25);
            threshold           = 0.05;
            [lct, pct, ~, ~, ~] = fx_track_byphase_2pi_increment(Idx_resin, Ascan_as, threshold);
            % drop the Nan in 'I_phase_track'
            I_phase_track       = lct(2:end-1);
            I_phase_track       = I_phase_track(~isnan(I_phase_track));
            P_phase_track       = pct(2:end-1);
            P_phase_track       = P_phase_track(~isnan(P_phase_track));
            % plot the inam and the record signal
            cf       = figure('Name', strcat('plytrack', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca       = subplot(2, 1, 1);
            plot(ca, (1:length(Ascan_interp))/obj.fs*1e6, Ascan_interp, 'linewidth', 2);
            hold on;
            plot(ca, (1:length(Ascan_as))/obj.fs*1e6, abs(Ascan_as), 'linewidth', 2);
            % use extracted location of rear surface here.
            plot(ca, [front_I_interp lct(end)]/obj.fs*1e6, abs(Ascan_as(round([front_I_interp lct(end)]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            hold on;
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            legend({'Recorded', 'Inst. Amp.'});
            legend('Location', 'best');
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % plot the inph, and the ply track of - pi / 2
            ca = subplot(2, 1, 2);
            plot(ca, (1:length(Ascan_as))/obj.fs*1e6, angle(Ascan_as), 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            % use extracted location of rear surface here.
            plot(ca, [front_I_interp lct(end)]/obj.fs*1e6, angle(Ascan_as(round([front_I_interp lct(end)]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            hold on;
            plot(ca, I_phase_track/obj.fs*1e6, P_phase_track, 'gd', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            legend({'Inst. phase', 'Derived surface', 'Derived inter-ply'});
            legend('Location', 'northwest');
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function demo_removefront_filter_plytrack(obj, x, y, f0, sigma)
            % demonstrate the A scan and the ply track after surface echo
            % removing and logGabor filter
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % n_interface: the number of the interfacec we want to track
            % f0: the array of the center frequencies of the filter
            % sigma: affects the bandwidth of the filter.
            Ascan = squeeze(obj.img(x, y, :));
            front_I_interp = obj.front_I(x, y) ;
            rear_I_interp = obj.rear_I(x, y) ;
            % find the origin inph_unw(1), used as the base line of front echo
            Ascan_as = hilbert(Ascan);
            A_inph =angle(Ascan_as(round(front_I_interp): round(rear_I_interp)));
            A_inph_unw = unwrap(A_inph);
            phi_0 = A_inph_unw(1);
            % judge who is longer, and process afterwards.
            lAs = length(Ascan);
            lref = length(obj.refer_Ascan_aligned);
            if lAs >= lref
                Ascan(1:lref) = Ascan(1:lref) - obj.refer_Ascan_aligned;
            else
                Ascan = Ascan - obj.refer_Ascan_aligned(1: lAs);
            end
            % define a logGabor filter and calcualte the analytical signal
            % calculate the AS in frequency domain
            Fo = fft(Ascan);
            % transverse Fo for dot multiply
            Fo = Fo';
            fx_shift = fftshift(Fo);
            n = length(Ascan);
            fshift = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            Fa = (1 + sign(fshift)) .* fx_shift;
            % add the logGabor filter here
            filter =  fx_1dLogGabor_filter(fshift, f0, sigma);
            Fa_filter = Fa .* filter;
            % ifft to obtain the AS
            Fa_ishift = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as = conj(flip(Ascan_as));
            % plot the inam and the record signal
            cf = figure('Name', strcat('plytrack_removingfront_filter', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(2, 1, 1);
            plot(ca, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, abs(Ascan_as), 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amp. (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Recorded', 'Inst. Amp.'})
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % track the - pi / 2 by interplation on the unwraped data
            A_inph = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            A_inph_unw = unwrap(A_inph);
            phase_seq = phi_0 + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
            I_phase_track = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                phase_seq) + front_I_interp;
            % plot the inph, and the ply track of - pi / 2
            ca = subplot(2, 1, 2);
            plot(ca, angle(Ascan_as), 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, I_phase_track, angle(Ascan_as(round(I_phase_track))) , 'bv', 'linewidth', 2);
            hold on;
            plot(ca, [front_I_interp rear_I_interp], angle(Ascan_as(round([front_I_interp rear_I_interp]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} ', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function show_track_interply(obj, xslice, yslice, zslice)
            % demonstrate the 3d instantaneous phase
            % show inph_ex
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            inph_ex = obj.mask_plytrack;
            y = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure('Name', ['3dtrack_interply_' num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            h = slice(ax, X, Y, Z, inph_ex, xslice, yslice, zslice);
            set(h, 'EdgeColor', 'none');
            colormap gray;
            %             h = colorbar;
            %             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. phase(rad.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            %
            %             [X, Y] = meshgrid(x, y);
            %             surf(ax, X, Y, obj.front_I / obj.fs * 1e6, 'FaceAlpha',0.5, 'FaceColor', 'red');
            %             hold on;
            %             surf(ax, X, Y, obj.rear_I / obj.fs * 1e6, 'FaceAlpha',0.5, 'FaceColor', 'magenta');
            %             hold on;
        end
        
        function [X, Y, profile_layer] = show_oneinterply(obj, layer, filtername, v, win_x, win_y, TOF)
            % demonstrate the profile of one inter-ply in 3d
            % layer: 'layer'th layer
            % filtername: 'nofilter' or 'logGabr', for naming the figure.
            % v: the average velocity of the laminate
            % win_x, win_y: the window to cut the 2d area
            % TOF: TOF of one ply, unit: sampling points
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            inph_ex = obj.mask_plytrack;
            if ~isempty(win_x)
                inph_ex = inph_ex(win_x, :, :);
            end
            if ~isempty(win_y)
                inph_ex = inph_ex(:, win_y, :);
            end
            profile_layer = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    k = find(inph_ex(i, j, :), layer);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                    if length(k) >= layer
                        profile_layer(i, j) = (k(end) - obj.front_I(i, j)) / obj.fs * 1e6;
                    else % set the last inter-ply as NaN otherwise
                        profile_layer(i, j) = NaN;
                    end
                end
            end
            figure('Name', ['one_interply_profile' , filtername, '_', num2str(layer)]);
            ax            = subplot(1, 1, 1);
            X             = (1: size(profile_layer, 1)) / obj.fx * 1e3;
            Y             = (1: size(profile_layer, 2)) / obj.fy * 1e3;
            % convert to distance mm from TOF us
            profile_layer = profile_layer * v / 2 / 1e3;
            imagesc(ax, X, Y, profile_layer);
            hold on;
            h = colorbar;
            % keep the color axis in one ply
            if TOF~=0
                ply_front     = layer*TOF-TOF/2;
                ply_back      = layer*TOF+TOF/2;
                caxis([ply_front/obj.fs*1e3*v/2, ply_back/obj.fs*1e3*v/2]); % for comparision, the 'caxis' needs to be modified
            end
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Depth (mm)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function show_B_scan_interply(obj, B_type, index, win, PropertyName)
            % demonstrate the B_scan of the interply track
            % show inph_ex
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % win: the win of the B-scan image
            % PropertyName: 'img_hil' or 'img_hil_filter' ... used for the inam background.
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            inph_ex     = logical(obj.mask_plytrack);
            x           = (1:(win(end)-win(1)+1))/obj.fx*1e3;
            y           = (1:(win(end)-win(1)+1))/obj.fx*1e3;
            inam_visual = abs(obj.(PropertyName));
            inph_visual = angle(obj.(PropertyName));
            z           = (1:size(inam_visual, 3))/obj.fs*1e6;
            cf          = figure('Name', ['B_scan_inam_' PropertyName, '_', B_type, '_', num2str(index)]);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            % 3D Gaussian filter
            sigma1          = [3 3 3];
            inph_cos        = cos(inph_visual);
%             inph_sin        = sin(inph_visual);
            inph_cos_filter = imgaussfilt3(inph_cos, sigma1, 'padding', 'replicate');
%             inph_sin_filter = imgaussfilt3(inph_sin, sigma1, 'padding', 'replicate');
            inph_visual     = acos(inph_cos_filter);
            if (B_type == 'x')
                % B scan image
                B_scan      = squeeze(inam_visual(win, index, :));
                B_scan_inph = squeeze(inph_visual(win, index, :));
                % front line
                front_line  = squeeze(obj.front_I(win, index));
                % rear line
                rear_line   = squeeze(obj.rear_I(win, index));
                % interply track
                if ~isempty(inph_ex)
                    inph_ex     = squeeze(inph_ex(win, index, :));
                end
            elseif (B_type == 'y')
                % B scan image
                B_scan      = squeeze(inam_visual(index, win, :));
                B_scan_inph = squeeze(inph_visual(index, win, :));
                % front line
                front_line  = squeeze(obj.front_I(index, win));
                % rear line
                rear_line   = squeeze(obj.rear_I(index, win));
                % interply track
                if ~isempty(inph_ex)
                    inph_ex     = squeeze(inph_ex(index, win, :));
                end
            end
%             % filter 
%             for i = 1:size(B_scan, 1)
%                 B_scan(i, :) = bandpass(B_scan(i, :), [10e6 20e6], obj.fs);
%             end
            front_line = front_line/obj.fs*1e6;
            rear_line  = rear_line/obj.fs*1e6;
            % choose the 'x' axis
            if (B_type == 'y')
                x = y;
            end
            imagesc(x, z, abs(B_scan'));
            hold on;
            % figure setup
            colormap(gray);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16}Inst. amp. (arb)');
%             cl   = caxis;
%             caxis(cl);
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            % draw the line indicating the front surface
            hold on;
            plot(x, front_line, 'r-', 'linewidth', 1);
            % draw the line indicating the rear surface
            hold on;
            plot(x, rear_line, 'magenta-', 'linewidth', 1);
            [row_inph, col_inph] = find(inph_ex);
            row_inph             = row_inph  / obj.fx * 1e3;
            col_inph             = col_inph / obj.fs * 1e6;
            hold on;
            scatter(row_inph, col_inph, 5, 'o', 'cyano', 'filled');
            legend({'Front surface', 'Rear surface', 'Interply track'});
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            % inph 
            cf = figure('Name', ['B_scan_inph_', PropertyName, '_', B_type, '_', num2str(index)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            imagesc(x, z, B_scan_inph');
            hold on;
            % figure setup
            colormap(gray);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16}Inst. phase. (rad.)');
            cl   = caxis;
            caxis(cl);
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            % draw the line indicating the front surface
            hold on;
            plot(x, front_line, 'r-', 'linewidth', 1);
            % draw the line indicating the rear surface
            hold on;
            plot(x, rear_line, 'magenta-', 'linewidth', 1);
            legend({'Front surface', 'Rear surface', 'Interply track'});
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function [TOF, x, z, inph_ex, Bscan_background, front_line, rear_line]= track_show_B_scan_interply(obj, B_type, index, win, PropertyName, threshold, nol, std_mean_threshold)
            % demonstrate the B_scan of the interply track
            % show inph_ex
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % win: the win of the B-scan image
            % PropertyName: 'img_hil' or 'img_hil_filter' ... used for the inam background.
            % nol: number of layers    
            % std_mean_threshold: the std and mean to judge the trackable
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            inph_visual = obj.(PropertyName);
            x           = (1:(win(end)-win(1)+1)) / obj.fx * 1e3;
            y           = (1:(win(end)-win(1)+1)) / obj.fy * 1e3;
            z           = (1: size(inph_visual, 3) ) / obj.fs * 1e6;   
            inph_ex     = zeros(length(x), length(z));
            thickness   = NaN(length(x), nol-1);
            % debug
            B_scan_true = zeros(length(x), length(z));
            % %             mask to select the interply track;
            %             for i = 1:size(inph_ex, 1)
            %                 for j = 1:size(inph_ex, 2)
            %                     indx(i, j, :) = find(inph_ex(i, j, :));
            %                 end
            %             end
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure('Name', ['B_scan_interplytrack_' PropertyName, '_', B_type, '_', num2str(index)]);
            %             ax = subplot(1, 1, 1);
            errors_pos_repeat  = NaN(length(x), nol);
            if (B_type == 'x')
                % B scan image
                B_scan     = squeeze(inph_visual(win, index, :));
                % front line
                front_line = squeeze(obj.front_I(win, index));
                % rear line
                rear_line  = squeeze(obj.rear_I(win, index));
            elseif (B_type == 'y')
                % B scan image
                B_scan     = squeeze(inph_visual(index, win, :));
                % front line
                front_line = squeeze(obj.front_I(index, win));
                % rear line
                rear_line  = squeeze(obj.rear_I(index, win));
            end
            % save and return the TOF of front- and back-walls echoes
            TOF = front_line - rear_line;
            % track the interface by phase
            for idx = 1:size(B_scan, 1)
                if rear_line(idx)<=0 || rear_line(idx)>size(B_scan, 2) ...
                        || front_line(idx)<=0 || front_line(idx)>rear_line(idx)
                    continue;
                end
                A_ori     = squeeze(B_scan(idx, :));
                Idx_resin = linspace(front_line(idx), rear_line(idx), 25);
                % track the -pi/2 by phase 
                power_noise            = 1;
                if strcmp(PropertyName, 'img_hil') || strcmp(PropertyName, 'img_hil_filter')
                    [lct, ~, err_p, ~, ~] = fx_track_byphase_step(Idx_resin, real(A_ori), obj.fs, threshold, power_noise);
                else % track by peaks
                    [lct, ~, err_p, ~]    = fx_track_bypeaks(Idx_resin, A_ori, threshold, power_noise);
                end
                % track the -pi/2 by distance
                %                 [lct, ~, ~]       = fx_track_bycluster(A_ori, threshold, 23, Idx_resin(2:end-1));
                %                 [lct, ~, ~, errors_pos]     = fx_track_bydistance(A_ori, threshold, Idx_resin);
                % check nan
                for idx_z = 2:length(lct)-1
                    if ~isnan(lct(idx_z)) % less than the rear surface
                        inph_ex(idx, round(lct(idx_z))) = 1;
                    end
                end
                errors_pos_repeat(idx, :) = err_p;
                % thickness calculation
                thickness(idx, :)         = -diff(lct);
                % update the back rear
                if isnan(lct(end))
                    rear_line(idx) = find(inph_ex(idx, :), 1, 'last'); % find last non zero value
                else
                    rear_line(idx) = lct(end); 
                end
                % debug
                B_scan_true(idx, round(Idx_resin(2:end-1))) = 1;
%                 front_line(idx)                    = lct(1);
%                 rear_line(idx)                     = lct(end);
            end
            front_line           = front_line / obj.fs  * 1e6;
            rear_line            = rear_line / obj.fs * 1e6;
            [row_inph, col_inph] = find(inph_ex);
            row_inph             = row_inph  / obj.fx * 1e3;
            col_inph             = col_inph / obj.fs * 1e6;
            % choose the 'x' axis
            if (B_type == 'y')
                x = y;
            end
            ca       = subplot(1, 1, 1);
            Bscan_background = abs(B_scan');
            imagesc(ca, x, z, Bscan_background);
            hold on;
            % figure setup
            colormap(gray);
            h = colorbar('eastoutside');
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16}Inst. amp. (arb)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            % draw the line indicating the front surface
            hold on;
            plot(ca, x, front_line, 'r-', 'linewidth', 2);
            % draw the line indicating the rear surface 
            hold on;
            scatter(ca, row_inph, col_inph, 8, 'cyano', 'd', 'filled');
            hold on;
            plot(ca, x, rear_line, 'magenta-', 'linewidth', 2);
            hold on;
            % debug
            [row_inph, col_inph] = find(B_scan_true);
            row_inph             = row_inph  / obj.fx * 1e3;
            col_inph             = col_inph / obj.fs * 1e6;
            hold on;
            scatter(ca, row_inph, col_inph, 3, [0.8,0.8,0.8], 'o', 'filled');
            hold on;
            legend({'Front surface', 'Estimated locations of interplies', 'Estimated rear surface', 'uniformly-spaced positions of interplies'});
            set(ca, 'fontsize', 16);
            set(ca, 'linewidth', 2);
            set(cf, 'Position', [-200, -200, 600, 600], 'color', 'white');
            %
            cf = figure('Name', ['Est_thickness' PropertyName, '_', B_type, '_', num2str(index)]);
            set(cf, 'Position', [-200, -200, 500, 600], 'color', 'white');
            %
            %
            ca2            = subplot(1, 1, 1);
            % plot the thickness by mean and var
            n              = nol - 1; 
            total_thick    = 5.52; % modify here
            thickness      = thickness / mean(TOF) * total_thick;      
            mean_thickness = squeeze(mean(thickness, 1, 'omitnan'));
            std_thickness  = squeeze(std(thickness, 1, 'omitnan'));
            untrackable    = mean_thickness>total_thick/n*(1+std_mean_threshold(1)) | ...
                mean_thickness<total_thick/n*(1-std_mean_threshold(1)) | ...
                std_thickness>total_thick/n*std_mean_threshold(2);
            h = fx_shaded_errorbars(ca2, untrackable, (1:n), mean_thickness, std_thickness);
            hold on;
            yl = yline(total_thick/n, 'k--', '', 'LineWidth',2);
            yl.LabelHorizontalAlignment = 'right';
            grid on;
            % legend
            legend({'Estimated thickness', 'Uniformly-spaced thickness'}, 'location', 'best');
            %             xticks(1:2:24);
            xlabel('\fontname {times new roman} layer i', 'fontsize', 16);
            ylabel('\fontname {times new roman} Thickness (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca2, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca2, 'linewidth', 2);
%             ylim([0.15 0.3]); % modify here
        end

        function show_Cscan_fitplytrack(obj, layer, ratio, PropertyName)
            % demonstrate the profile of one inter-ply in 3d
            % units: mm, mm, us
            % layer: the number of the layer
            % ratio: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % ***
            % the ratio should be 0 ~ 1
            if ratio < 0 || ratio > 1
                msg = "The 'raio' should be between 0 and 1";
                error(msg);
            end
            inph_ex = obj.mask_plytrack;
            inph = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            inam = abs(obj.(PropertyName));
            profile_layer = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_inph = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_inam = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    inph_ex(i, j, [round(obj.front_I(i, j)) round(obj.rear_I(i, j))]) = 1; % add the front_surface and the rear_surface
                    k = find(inph_ex(i, j, :), layer + 1);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                    % extract two inter-plies besides the ply
                    if length(k) > layer
                        k1 = k(layer);
                        k2 = k(layer + 1);
                        % use the ratio to calculate the position in the ply
                        z_index = round(k1 + (k2 - k1) * ratio);
                        C_scan_inph(i, j) = inph(i, j, z_index);
                        C_scan_inam(i, j) = inam(i, j, z_index);
                        profile_layer(i, j) = z_index / obj.fs * 1e6;
                    else % set the last inter-ply as NaN otherwise
                        profile_layer(i, j) = NaN;
                    end
                end
            end
            figure('Name', ['Cscan_inam_fitplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio)]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan_inam, 2)) / obj.fy * 1e3;
            imagesc(ax, X,Y, C_scan_inam)
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % plot the profile of the layer
            figure('Name', ['profile_fitplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio)]);
            ax = subplot(1, 1, 1);
            X = (1: size(profile_layer, 1)) / obj.fx * 1e3;
            Y = (1: size(profile_layer, 2)) / obj.fy * 1e3;
            imagesc(ax, X,Y, profile_layer)
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            %
            %             figure('Name', ['Cscan_inph_fitplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio)]);
            %             ax = subplot(1, 1, 1);
            %             X = (1: size(C_scan_inph, 1)) / obj.fx * 1e3;
            %             Y = (1: size(C_scan_inph, 2)) / obj.fy * 1e3;
            %             imagesc(ax, X,Y, C_scan_inph)
            %             hold on;
            %             h = colorbar;
            %             %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            %             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            %             set(ax, 'fontsize', 16);
            %             set(ax, 'linewidth', 1.5);
            %             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            %             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(gca, 'linewidth', 2);
        end
        
        function compare_Cscan_fitplytrack(obj, layer, ratio, PropertyName)
            % demonstrate the profile of one inter-ply in 3d
            % units: mm, mm, us
            % layer: the number of the layer
            % ratio: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % ***
            % the ratio should be 0 ~ 1
            if ratio < 0 || ratio > 1
                msg = "The 'raio' should be between 0 and 1";
                error(msg);
            end
            inph_ex          = obj.mask_plytrack;
            inph             = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            inam             = abs(obj.(PropertyName));
            profile_layer    = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_inph      = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_inam_temp = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    inph_ex(i, j, [round(obj.front_I(i, j)) round(obj.rear_I(i, j))]) = 1; % add the front_surface and the rear_surface
                    k = find(inph_ex(i, j, :), layer + 1);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                    % extract two inter-plies besides the ply
                    if length(k) > layer
                        k1                     = k(layer);
                        k2                     = k(layer + 1);
                        % use the ratio to calculate the position in the ply
                        z_index                = round(k1 + (k2 - k1) * ratio);
                        C_scan_inph(i, j)      = inph(i, j, z_index);
                        C_scan_inam_temp(i, j) = inam(i, j, z_index);
                        profile_layer(i, j)    = z_index / obj.fs * 1e6;
                    else % set the last inter-ply as NaN otherwise
                        profile_layer(i, j)    = NaN;
                    end
                end
            end
            figure('Name', ['Cscan_inam_fitplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio)]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_temp, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan_inam_temp, 2)) / obj.fy * 1e3;
            imagesc(ax, X,Y, C_scan_inam_temp)
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % plot the parallel C-scan of the layer
            figure('Name', ['Cscan_inam_' , PropertyName, '_', num2str(z_index)]);
            ax = subplot(1, 1, 1);
            z_index = round(mean(profile_layer * obj.fs / 1e6, 'all', 'omitnan'));
            C_scan_inam_2 = inam(:, :, z_index);
            C_scan_inam_2(isnan(C_scan_inam_temp)) = NaN;
            imagesc(ax, X,Y, C_scan_inam_2);
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function [mean_thickness, std_thickness] = thickness_estimation_nolimit(obj, win_x, win_y, total_thick, figname)
            % calculate the thickness of each ply
            % demonstrate the thickness by mean and variance
            % do no assumpt the number of layers.
            % win_x, win_y: the window to cut the 2d area
            % total_thick: the thickness of the sample
            % std_mean_threshold: threshold to determine the trackable
            % figname: name of the figure
            TOF            = mean(obj.rear_I - obj.front_I, 'all'); % total samples
            thickness      = obj.est_thickness;
            thickness      = thickness(win_x, win_y, :)/TOF *total_thick;
            mean_thickness = squeeze(mean(thickness, [1, 2], 'omitnan'));
            std_thickness  = squeeze(std(thickness, 0, [1, 2], 'omitnan'));
            % plot the thickness by mean and var
            figure('Name', ['thickness_estimation_' , figname, '_', ...
                num2str(total_thick * 100)]);
            ca2 = subplot(1, 1, 1);
             % plot the thickness by mean and var
            untrackable    = ~find(mean_thickness);
            h = fx_shaded_errorbars(ca2, untrackable, 1:length(mean_thickness), mean_thickness, std_thickness);
            hold on;
            grid on;
            % legend
            legend({'Estimated thickness', 'Uniformly-spaced thickness'}, 'location', 'best');
            %             xticks(1:2:24);
            xlabel('\fontname {times new roman} Ply i', 'fontsize', 16);
            ylabel('\fontname {times new roman} Thickness (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca2, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca2, 'linewidth', 2);
        end
        
        function [mean_thickness, std_thickness] = thickness_estimation(obj, total_thick, n, win_x, win_y, std_mean_threshold, figname)
            % calculate the thickness of each ply
            % demonstrate the thickness by mean and variance
            % total_thick: the total thickness of the composite unit: mm
            % n: the number of the layers to be found
            % win_x, win_y: the window to cut the 2d area
            % std_mean_threshold: threshold to determine the trackable
            % figname: name of the figure
            TOF            = mean(obj.rear_I - obj.front_I, 'all'); % total samples
            thickness      = obj.est_thickness;
            thickness      = thickness(win_x, win_y, :)/TOF *total_thick;
            mean_thickness = squeeze(mean(thickness, [1, 2], 'omitnan'));
            std_thickness  = squeeze(std(thickness, 0, [1, 2], 'omitnan'));
            % plot the thickness by mean and var
            figure('Name', ['thickness_estimation_' , figname, '_', ...
                num2str(total_thick * 100), '_', num2str(n)]);
            ca2 = subplot(1, 1, 1);
             % plot the thickness by mean and var
            total_thick    = 5.5; % modify here
            untrackable    = mean_thickness>total_thick/n*(1+std_mean_threshold(1)) | ...
                mean_thickness<total_thick/n*(1-std_mean_threshold(1)) | ...
                std_thickness>total_thick/n*std_mean_threshold(2);
            h = fx_shaded_errorbars(ca2, untrackable, (1:n), mean_thickness, std_thickness);
            hold on;
            yl = yline(total_thick/n, 'k--', '', 'LineWidth',2);
            yl.LabelHorizontalAlignment = 'right';
            grid on;
            % legend
            legend({'Estimated thickness', 'Uniformly-spaced thickness'}, 'location', 'best');
            %             xticks(1:2:24);
            xlabel('\fontname {times new roman} Ply i', 'fontsize', 16);
            ylabel('\fontname {times new roman} Thickness (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca2, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca2, 'linewidth', 2);
        end
        
        function [mean_thickness, var_thickness] = thickness_estimation_v1(obj, total_thick, n, win_x, win_y, figname)
            % calculate the thickness of each ply
            % demonstrate the thickness by mean and variance
            % total_thick: the total thickness of the composite unit: mm
            % n: the number of the layers to be found
            % win_x, win_y: the window to cut the 2d area
            % figname: name of the figure
            inph_ex   = obj.mask_plytrack;
            inph_ex   = inph_ex(win_x, win_y, :);
            thickness = NaN(size(inph_ex, 1), size(inph_ex, 2), n);
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    ls =  obj.rear_I(i, j) -  obj.front_I(i, j); % total samples
                    k  = find(inph_ex(i, j, :), n-1);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                    % 1st
                    if k(1) > obj.front_I(i, j)
                        thickness(i, j, 1) = (k(1) - obj.front_I(i, j)) / ls * total_thick ;
                    end
                    thickness(i, j, 2:length(k)) = diff(k) / ls * total_thick ;
                    % last
                    if obj.rear_I(i, j) > k(end)
                        thickness(i, j, n) = (obj.rear_I(i, j) - k(end)) / ls * total_thick ;
                    end
                end
            end
            mean_thickness = squeeze(mean(thickness, [1, 2], 'omitnan'));
            var_thickness  = squeeze(std(thickness, 0, [1, 2], 'omitnan'));
            % plot the thickness by mean and var
            figure('Name', ['thickness_estimation_' , figname, '_', ...
                num2str(total_thick), '_', num2str(n)]);
            ax = subplot(1, 1, 1);
            errorbar(1:n, mean_thickness*1e3, var_thickness*1e3, ...
                's', 'MarkerSize', 5, 'MarkerEdgeColor', 'red', ...
                'MarkerFaceColor','red', 'CapSize', 5, 'LineWidth', 2, 'LineStyle', 'none')
            hold on;
            yl = yline(total_thick/n, 'k--', ...
                'uniformly-spaced thickness', 'LineWidth',2);
            yl.LabelHorizontalAlignment = 'right';
            grid on;
            xticks(1:2:24);
            xlabel('\fontname {times new roman} CFRP ply {\itp}', 'fontsize', 16);
            ylabel('\fontname {times new roman} Thickness (\mum)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            box on;
        end
        
        function obj = error_pos_calculate(obj, win_x, win_y, layers_num, PropertyName)
            % calcalate the position errors
            % win_x, win_y: the window to cut the 2d area 
            % layers_num: the number of the layers to track
            % PropertyName: the name of the data property
            % ************
            mask_plytrack_temp = obj.mask_plytrack;
            front_I_temp       = obj.front_I;
            rear_I_temp        = obj.rear_I;
            errors_pos_repeat  = NaN(size(mask_plytrack_temp,1), ...
                size(mask_plytrack_temp,2), layers_num);
            for idx = win_x 
                for idy = win_y
                    if rear_I_temp(idx, idy)<=0 || rear_I_temp(idx, idy)>size(mask_plytrack_temp, 3) ...
                            || front_I_temp(idx, idy)<=0 || front_I_temp(idx, idy)>rear_I_temp(idx, idy)
                        continue;
                    end
                    Idx_resin            = linspace(front_I_temp(idx, idy), ...
                        rear_I_temp(idx, idy), layers_num);  % !!!!! 25 is the number of layers incl. front and back surfaces
                    A_mask               = squeeze(mask_plytrack_temp(idx, idy, :));
                    lct                  = find(A_mask);
                    lct_incl_suf         = [front_I_temp(idx, idy) lct' rear_I_temp(idx, idy)];
                    if length(lct_incl_suf)==length(Idx_resin)
                        % find the nearest one to calculate the errors
                        for idx_l = 1:length(Idx_resin)
                            [~, I_err]                         = min(abs(lct_incl_suf - Idx_resin(idx_l)));
                            errors_pos_repeat(idx, idy, idx_l) = (lct_incl_suf(I_err) - Idx_resin(idx_l))/obj.fs*1e6;
                        end
%                         errors_pos_repeat(idx, idy, :) = (lct_incl_suf-Idx_resin) / obj.fs * 1e6; % unit: us;
                    else
                        % find the nearest one to calculate the errors
                        for idx_l = 1:length(Idx_resin)
                            [~, I_err]                         = min(abs(lct_incl_suf-Idx_resin(idx_l)));
                            errors_pos_repeat(idx, idy, idx_l) = (lct_incl_suf(I_err)-Idx_resin(idx_l))/obj.fs*1e6;
                        end
%                         errors_pos_repeat(idx, idy, 1:length(lct_incl_suf)) = ...
%                             (lct_incl_suf-Idx_resin(1:length(lct_incl_suf))) / obj.fs * 1e6; % unit: us;
%                         errors_pos_repeat(idx, idy, 1:length(lct_incl_suf)-1) = (lct_incl_suf(2:end)...
%                             - lct_incl_suf(1:end-1) - Idx_oneply) / obj.fs * 1e6; % unit: us;
                    end
                end
                disp([num2str(idx), '/' , num2str(max(win_x))]);
            end    
            obj.errors_pos    = errors_pos_repeat;
            % error bar
            error_mean = squeeze(mean(abs(errors_pos_repeat),[1 2], 'omitnan'));
            error_std  = squeeze(std(errors_pos_repeat, 0, [1 2], 'omitnan'));
            cf = figure('Name', ['error_bar' '_' PropertyName]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(1, 1, 1);
            errorbar((1:layers_num)/obj.fs*1e6, error_mean*100, error_std*100,'-s',...
                'MarkerSize', 10, 'Linewidth', 3,'MarkerEdgeColor','red','MarkerFaceColor','red')
            xlabel('\fontname{times new roman} Time (µm)', 'fontsize', 16);
            ylabel('\fontname{times new roman} Errors (mean and std)', 'fontsize', 16);
            ytickformat('percentage');
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            grid on;
        end
        
        function makeMovie_InterplyTrack(obj, B_type, win, property_name)
            % make a movie of the B_scan of the interply track
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % win: the win of the B-scan image
            % property_name: name of the dataset
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            dataset = abs(obj.(property_name));
            front_p = obj.front_I;
            bakc_p  = obj.rear_I;
            x       = (1:(win(end)-win(1)+1))/obj.fx * 1e3;
            y       = (1:(win(end)-win(1)+1))/obj.fy * 1e3;
            z       = (1:size(dataset, 3))/obj.fs * 1e3 * 3000/2;
            cf      = figure('Name', ['B_scan_angles_', B_type, '_', num2str(win(end))]);
            set(cf, 'Position', [0, 0, length(win)*2, 1200], 'color', 'white');
            % make vedio
            axis tight manual;
            title(['y =' num2str(0/obj.fy * 1e3) ' mm']);
            xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} z (mm)', 'fontsize', 16);
            v = VideoWriter([property_name, '.avi']);
            open(v);
            %
            if (B_type == 'x')
                for index = 1:size(dataset,1)
                    % B scan image
                    % side view
                    B_scan = squeeze(dataset(index, win, :));
                    figure(cf);
                    ca1 = subplot(2, 1, 1);
                    % imagesc(x, z, B_scan');
                    pcolor(x, z, B_scan');
                    shading flat;
                    set(ca1, 'ydir', 'reverse');
                    hold on;
                    title(['y =' num2str(index/obj.fy * 1e3) ' mm']);
                    xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
                    ylabel('\fontname {times new roman} z (mm)', 'fontsize', 16);
                    mask  = (obj.row==index) & obj.col>=win(1) & obj.col<=win(end);
                    x_row = obj.col(mask) / obj.fx * 1e3;
                    x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
                    h1 = scatter(ca1, x_row, x_dep, ...
                        3,  'cyan', 'filled', ...
                        'DisplayName','Interply track');
                    hold on;
                    h2 = scatter(x, front_p(index, win)/ obj.fs * 1e3 * 3000/2, ...
                        3, 'red', 'filled', ...
                        'DisplayName','Front surface');
                    hold on;
                    h3 = scatter(x, bakc_p(index, win)/ obj.fs * 1e3 * 3000/2, ...
                        3, 'magenta', 'filled', ...
                        'DisplayName','Rear surface');
                    hold on;
                    % figure setup
                    colormap(gray);
                    caxis(ca1, [0 0.6]);
                    h = colorbar;
                    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
                    set(ca1, 'Fontname', 'times new Roman');
                    set(ca1, 'Fontsize', 16);
                    set(ca1, 'Linewidth', 1.5);
                    % top view
                    figure(cf);
                    ca2     = subplot(2, 1, 2);
                    % imagesc(x, z, B_scan');
                    z_index     = round(mean(front_p, 'all', 'omitnan') ...
                        + mean(bakc_p-front_p, 'all', 'omitnan') / size(dataset,1) * index);
                    B_scan = squeeze(dataset(1:end, 1:end, z_index));
                    pcolor(x, y, B_scan);
                    shading flat;
                    hold on;
                    title(['z =' num2str(z_index/obj.fs * 1e3 * 3000/2) ' mm']);
                    xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
                    ylabel('\fontname {times new roman} y (mm)', 'fontsize', 16);
                    hold on;
                    % figure setup
                    colormap(ca2, jet);
                    caxis(ca2, [0 0.6]);
                    h = colorbar;
                    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
                    set(ca2, 'Fontname', 'times new Roman');
                    set(ca2, 'Fontsize', 16);
                    set(ca2, 'Linewidth', 1.5);
                    % frame vedio operation
                    frame = getframe(cf);
                    writeVideo(v, frame);
                    clf(cf);
                end
            elseif (B_type == 'y')
            end
        end
        
        % ******************* reference signal **************
        function obj = read_refer(obj, filename, x, y)
            % read the reference signal, which is produced by one Ascan
            % filename: the pathname + filename of the .tdms file
            % x, y: the [x, y] idx to select the reference signal
            obj.reference = class_reference_signal(filename);
            obj.refer_Ascan_ori = squeeze(obj.reference.oimg_uint8(x, y, :)); % for reference signal , the index is normally (1, 1)
        end
        
        function show_reference_signal(obj)
            % display the Ascan, and its hilbert, freq. from the origin
            % ref. object
            % for reference signal , the index is normally (1, 1)
            obj.reference.show_hilbert_ref(1, 1);
        end
        
        function obj = cut_reference_signal(obj, win)
            % cut the reference signal if there are extra peaks in it.
            % the points outsides the window will be set 0
            % ********
            % the unit of 'win' is time second. Transfer to index
            idx1 = max(1, round( win(1) * obj.reference.fs));
            idx2 = round( win(2) * obj.reference.fs);
            obj.refer_Ascan = obj.refer_Ascan_ori(idx1: idx2);
            L = length(obj.refer_Ascan);
            n = 2^nextpow2(L);
            Y = fft(obj.refer_Ascan, n);
            f = obj.reference.fs .* (0:(n / 2)) / n;
            P = abs(Y / n);
            figure();
            subplot(2, 1, 1);
            plot(obj.refer_Ascan);
            subplot(2, 1, 2);
            plot(f,P(1:n/2+1));
        end
        
        function obj = align_refer_ascan(obj, x, y)
            % align the reference ascan and the recorded ascan
            % x, y: the x, y index to select the A scan
            record_signal = squeeze(obj.img(x, y, :));
            record_signal_hil = hilbert(record_signal);
            refer_Ascan_hil = hilbert(obj.refer_Ascan);
            t =  (1:length(record_signal)) / obj.reference.fs;
            t_ref = (1:length(obj.refer_Ascan)) / obj.reference.fs; % the time space of reference signal
            % display the figures
            cf = figure('Name', 'ref_record_signal');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            plot(t, record_signal, 'b-');
            hold on;
            plot(t, abs(record_signal_hil), 'b--');
            hold on;
            plot(t_ref, obj.refer_Ascan, 'r-');
            hold on;
            plot(t_ref, abs(refer_Ascan_hil), 'r--');
            xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
            ylabel('\fontname {times new roman}', 'fontsize', 16);
            title('real component', 'fontsize', 16, 'Fontname', 'times new Roman');
            hold on;
            % find the index and amp. of the max in inam
            [C_record, I_record] = max(abs(record_signal_hil));
            [C_ref, I_ref] = max(abs(refer_Ascan_hil));
            % shift the ref. signal by the indexes
            refer_Ascan_shift = circshift(obj.refer_Ascan, I_record - I_ref);
            refer_Ascan_shift_nor = refer_Ascan_shift * C_record / C_ref;
            % display the aligned ref. signal
            cf = figure('Name', 'alignref_record_signal');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            plot(t, record_signal, 'b-', 'linewidth', 2);
            hold on;
            plot(t_ref, refer_Ascan_shift_nor, 'r--', 'linewidth', 3);
            xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
            ylabel('\fontname {times new roman}', 'fontsize', 16);
            title('real component', 'fontsize', 16, 'Fontname', 'times new Roman');
            hold on;
            % save the aligned ref. signal in properties.
            obj.refer_Ascan_aligned = refer_Ascan_shift_nor;
        end
        
        function demo_removefront_plytrack(obj, x, y)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % interp_factor: the interpolation multiple
            % n_interface: the number of the interfacec we want to track
            Ascan = squeeze(obj.img(x, y, :));
            front_I_interp = obj.front_I(x, y) ;
            rear_I_interp = obj.rear_I(x, y) ;
            % find the origin inph_unw(1), used as the base line of front echo
            Ascan_as = hilbert(Ascan);
            A_inph =angle(Ascan_as(round(front_I_interp): round(rear_I_interp)));
            A_inph_unw = unwrap(A_inph);
            phi_0 = A_inph_unw(1);
            % judge who is longer, and process afterwards.
            lAs = length(Ascan);
            lref = length(obj.refer_Ascan_aligned);
            if lAs >= lref
                Ascan(1:lref) = Ascan(1:lref) - obj.refer_Ascan_aligned;
            else
                Ascan = Ascan - obj.refer_Ascan_aligned(1: lAs);
            end
            % plot the inam and the record signal
            Ascan_as = hilbert(Ascan);
            cf = figure('Name', ['removefront_plytrack' '_' num2str(x), '_', num2str(y)]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(2, 1, 1);
            plot(ca, Ascan, 'linewidth', 2);
            hold on;
            plot(ca, abs(Ascan_as), 'linewidth', 2)
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amp. (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Recorded', 'Inst. Amp.'})
            % track the - pi / 2 by interplation on the unwraped data
            A_inph = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            A_inph_unw = unwrap(A_inph);
            phase_seq = phi_0 + 3 / 2 * pi :  2 * pi : A_inph_unw(end);
            I_phase_track = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                phase_seq) + front_I_interp;
            % drop the Nan in 'I_phase_track'
            I_phase_track = I_phase_track(~isnan(I_phase_track));
            % plot the inph, and the ply track of - pi / 2
            ca = subplot(2, 1, 2);
            plot(ca, angle(Ascan_as), 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, I_phase_track, angle(Ascan_as(round(I_phase_track))) , 'bv', 'linewidth', 2);
            hold on;
            plot(ca, [front_I_interp rear_I_interp], angle(Ascan_as(round([front_I_interp rear_I_interp]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} ', 'fontsize', 16);
        end
        
        function r = calculate_SNR(obj, win_s, win_n)
            % obtain the SNR by calculating the energy in different window with the same length
            % win_s: the win to select the signal
            % win_n: the win to select the noise
            idx_s1 = round( win_s(1) * obj.reference.fs);
            idx_s2 = round( win_s(2) * obj.reference.fs);
            idx_n1 = round( win_n(1) * obj.reference.fs);
            idx_n2 = round( win_n(2) * obj.reference.fs);
            x      = obj.refer_Ascan(idx_s1: idx_s2);
            y      = obj.refer_Ascan(idx_n1: idx_n2);
            r      = snr(x, y);
            cf = figure('Name', ['reference_signal' '_' num2str(win_s(1)), '_', num2str(win_n(1))]);
            set(cf, 'Position', [0, 0, 1000, 600], 'color', 'white');
            t_space = (1:length(obj.refer_Ascan)) / obj.reference.fs * 1e6;
            ca = subplot(1, 1, 1);
            plot(ca, t_space, obj.refer_Ascan, 'LineWidth',2);
            hold on;
            rectangle('Position',[win_s(1)*1e6, min(obj.refer_Ascan), ...
                (win_s(2)-win_s(1))*1e6, max(obj.refer_Ascan)-min(obj.refer_Ascan)], ...
                'EdgeColor','g', 'LineWidth',2);
            hold on;
            text(win_s(1)*1e6, mean(obj.refer_Ascan)+ ...
                max(obj.refer_Ascan)/3-min(obj.refer_Ascan)/3 ...
                ,'Signal','Color','g','FontSize',16);
            hold on;
            rectangle('Position',[win_n(1)*1e6, min(obj.refer_Ascan), ...
                (win_n(2)-win_n(1))*1e6, max(obj.refer_Ascan)-min(obj.refer_Ascan)], ...
                'EdgeColor','r', 'LineWidth',2);
            hold on;
            text(win_n(1)*1e6, mean(obj.refer_Ascan)+ ...
                max(obj.refer_Ascan)/3-min(obj.refer_Ascan)/3 ...
                ,'Noise','Color','r','FontSize',16);
            xlabel('\fontname{times new roman} Time (\mus)', 'fontsize', 16);
            ylabel('\fontname{times new roman} Amp. (arb.)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
        end
        
        % ******************* ply track by peaks *********
        function obj = demo_barErrors_phase(obj, PropertyName, layers_num, threshold)
            % track the layers by peak 
            % display the error bar
            % PropertyName: the dataset to select, 'img', 'img_hil', 'img_deconv'
            % layers_num: the number of the layers to track
            % threshold:  threshold to find peaks
            img_ori             = obj.(PropertyName);
%             lct_repeat          = NaN(size(img_ori,1), size(img_ori,2), layers_num);
            errors_pos_repeat   = NaN(size(img_ori,1), size(img_ori,2), layers_num);
            for i = 1:size(img_ori,1)
                for j = size(img_ori,2)
                    ori_signal             = squeeze(img_ori(i, j, :));
                    Idx_resin              = linspace(obj.front_I(i, j), ...
                        obj.rear_I(i, j), layers_num);
                    % track the plies by peaks
                    [~, ~, errors_pos, ~]      = fx_track_byphase(Idx_resin, ori_signal, threshold);
                    %                 lct_repeat(i, :)        = lct;
                    errors_pos_repeat(i, j, :) = errors_pos;
                end
                disp(i);
            end
            % error bar
            error_mean = squeeze(mean(errors_pos_repeat,[1 2], 'omitnan'));
            error_std  = squeeze(std(errors_pos_repeat, 0, [1 2], 'omitnan'));
            cf = figure('Name', ['error_bar' '_' 'bypahse', '_', PropertyName]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(1, 1, 1);
            errorbar(1:layers_num, 100*error_mean, 100*error_std,'-s',...
                'MarkerSize', 10, 'Linewidth', 3,'MarkerEdgeColor','red','MarkerFaceColor','red')
            xlabel('\fontname{times new roman} Layer', 'fontsize', 16);
            ylabel('\fontname{times new roman} Errors (mean and std)', 'fontsize', 16);
            ytickformat('percentage');
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            grid on;
        end
        
        % ******************* deconvolution **************
        function obj = Wiener_deconv(obj, x, y, q_factor)
            % apply Wiener deconvolution on the A scan
            % x, y: the x, y index to select the A scan
            % q_factor: the Q factor for Wiener deconvolution
            record_signal = squeeze(obj.img(x, y, :));
            % move the maximum point to the center of the kernel
            kernel = obj.refer_Ascan;
            [~, I] = max(kernel);
            kernel = kernel(1: 2 * I);
            % deconvolve
            deconvolved = fx_wiener_deconvolution_1d(record_signal, kernel, q_factor);
            % plot origin
            t_space = (1:length(record_signal)) / obj.fs * 1e6;
            cf = figure;
            ca = subplot(2, 1, 1);
            plot (ca, t_space, record_signal, 'linewidth', 2);
            hold on;
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(2, 1, 2);
            plot(ca, t_space, deconvolved, 'linewidth', 2);
            hold on;
            ylabel(ca,{['\fontname {times new roman}\fontsize {16}'  'layers'] , '\fontname {times new roman}\fontsize {16} Amp. (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time (\mus)', 'fontsize', 16);
        end
        
        function obj = deconv_ARextrapolation(obj, B_win, B_type, index, q_factor_AR, fft_padding, bw, bwr, k, ker_wid, fig_subname)
            % apply the Wiener deconvolution combined with AR extrapolation on the original signals
            % Bwin: the indexes of the first and the last points in Bscan image
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % q_factor_AR: the Q_factor for wiener deconvolution combined with AR extrapolation;
            % fft_padding: the power of 2 for the length of the padded signal
            % bw: the bandwidht in unit dB
            % bwr: the the reference to define the bandwidht in unit dB
            % k: the order of AR regression
            % ker_wid: half of the width of the kernel, maximum centered
            % fig_subname: the subname for the fig
            % ****
            % move the maximum point to the center of the kernel
            kernel = obj.refer_Ascan;
            [~, I] = max(kernel);
            kernel = kernel(I - ker_wid: I + ker_wid);
            L = length(kernel);
            n = 2^nextpow2(L);
            Y = fft(kernel, n);
            f = obj.fs *(0:(n / 2)) / n;
            P = abs(Y(1:n/2+1) / n);
            if bwr~=0 % define bw by bwr
                [f1, f2] = fx_calculate_bw(f, P, bwr);
                bw = [f1, f2];
            end
            % create matrix
            [lx, ly, lz] = size(obj.img);
            x = (1: lx) / obj.fx * 1e3;
            y = (1: ly) / obj.fy * 1e3;
            z = (1: lz) / obj.fs * 1e6;
            if (B_type == 'x')
                % B scan image
                img_WienerDeconv_AR_temp = NaN(length(B_win), lz);
                for i = B_win
                    ori_signal = squeeze(obj.img(i, index, :));
                    % Wiener deconv. + ARregression
                    deconv = class_deconv(ori_signal, kernel, obj.fs);
                    wienerdeconv_ARextrap_s = deconv.wienerdeconv_ARextrap(q_factor_AR, bw, k, 'bg', fft_padding);
                    img_WienerDeconv_AR_temp(i - min(B_win) + 1, :) = wienerdeconv_ARextrap_s(1:lz);
                    disp([num2str(i - min(B_win) + 1) '/' num2str(lx)]);
                end
            elseif (B_type == 'y')
                img_WienerDeconv_AR_temp = NaN(length(B_win), lz);
                for i = B_win
                    ori_signal = squeeze(obj.img(index, i, :));
                    % Wiener deconv. + ARregression
                    deconv = class_deconv(ori_signal, kernel, obj.fs);
                    wienerdeconv_ARextrap_s = deconv.wienerdeconv_ARextrap(q_factor_AR, bw, k, 'bg', fft_padding);
                    img_WienerDeconv_AR_temp(i - min(B_win) + 1, :) = wienerdeconv_ARextrap_s(1:lz);
                    disp([num2str(i - min(B_win) + 1) '/' num2str(ly)]);
                end
            end
            obj.img_WienerDeconv_AR_B = img_WienerDeconv_AR_temp;
            % choose the 'x' axis
            if (B_type == 'y')
                x = y;
            end
            cf = figure('Name', ['Bscan_WienerDeconv_ARextrap' '_', B_type, '_', num2str(index), fig_subname]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            imagesc(x, z, img_WienerDeconv_AR_temp');
            % figure setup
            colormap(gray);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            caxis('auto');
            cl = caxis;
%             caxis(cl / 10);
        end
        
        function obj = apply_deconvolutions_Bscan(obj, B_win, B_type, index, q_factor, q_factor_AR, fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit)
            % apply the deconvolutions on the original signals
            % Bwin: the indexes of the first and the last points in Bscan image
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % Bwin: the indexes of the first and the last points in Bscan
            % q_factor: the Q_factor for wiener deconvolution, 1e-2 for  default;
            % q_factor_AR: the Q_factor for wiener deconvolution combined with AR extrapolation;
            % fft_padding: the power of 2 for the length of the padded signal
            % bw: the bandwidht in unit dB
            % bwr: the the reference to define the bandwidht in unit dB
            % k: the order of AR regression
            % lam: the regularization parameter for sparse deconvolution
            % ker_wid: half of the width of the kernel, maximum centered
            % Nit: the number of iteration for sparse deconvolution
            % DownRate: the parameter for downsampleing the kernel
            % fig_subname: the subname for the fig
            % ****
            % move the maximum point to the center of the kernel
            kernel = obj.refer_Ascan;
            [~, I] = max(kernel);
            kernel = kernel(I - ker_wid: I + ker_wid);
            % remove the direct component
            kernel = kernel - mean(kernel);
            L      = length(kernel);
            n      = 2^nextpow2(L);
            Y      = fft(kernel, n);
            f      = obj.fs *(0:(n / 2)) / n;
            P      = abs(Y(1:n/2+1) / n);
            % make a series of f_windows
            f_windows                    = NaN(length(bwr), 2);
            for idx_bwr = 1:length(bwr)
                bwr_single               = bwr(idx_bwr);
                [f1, f2]                 = fx_calculate_bw(f, P, bwr_single);
                f_windows(idx_bwr, :)    = [f1, f2];
            end
            % create matrix
            [lx, ly, lz] = size(obj.img);
            % B scan image
            img_WienerDeconv_temp    = NaN(length(B_win), lz);
            img_WienerDeconv_AR_temp = NaN(length(B_win), lz);
            img_SparseDeconv_temp    = NaN(length(B_win), lz);
            if (B_type == 'x')
                % B scan image
                for i = B_win
                    ori_signal                                      = squeeze(obj.img(i, index, :));
                    % remove the direct component
                    ori_signal                                      = ori_signal - mean(ori_signal);
                    % wiener deconvolve
                    wienerdeconv_s                                  = fx_wiener_deconvolution_1d(ori_signal', kernel', q_factor);
                    img_WienerDeconv_temp(i - min(B_win) + 1, :)    = abs(hilbert(wienerdeconv_s));
                    % Wiener deconv. + ARregression
                    deconv                                          = class_deconv(ori_signal, kernel, obj.fs);
                    wienerdeconv_ARextrap_s                         = deconv.wienerdeconv_ARextrap(q_factor_AR, f_windows, k, 'bg', fft_padding);
                    img_WienerDeconv_AR_temp(i - min(B_win) + 1, :) = abs(hilbert(wienerdeconv_ARextrap_s(1:lz)));
%                     % sparse deconvolution
%                     % downsample the kernel
%                     kernel_sparse = kernel(1:DownRate:end);
%                     L = length(kernel_sparse);
%                     M = length(ori_signal);
%                     N = M + 1 - L;               % M : length of observed signal
%                     % Create convolution matrix H
%                     H = sparse(M, N);
%                     e = ones(N, 1);
%                     for ii = 0:L-1
%                         H = H + spdiags(kernel_sparse(ii + 1) * e, -ii, M, N);            % H : convolution matrix (sparse)
%                     end
%                     %                     issparse(H)                                       % confirm that H is a sparse matrix
%                     % L1 norm penalty The penalty function is phi(x) = lam abs(x)
%                     phi_L1 = @(x) lam * abs(x);
%                     wfun_L1 = @(x) abs(x)/lam;
%                     %                     dphi_L1 = @(x) lam * sign(x);
%                     [x1, ~] = fx_deconv_MM(ori_signal, phi_L1, wfun_L1, H, Nit);
%                     img_SparseDeconv_temp(i - min(B_win) + 1, :) = [zeros(round(L / 2 - 1), 1); x1; zeros(round(L / 2) - 1, 1)];
                    disp([num2str(i - min(B_win) + 1) '/' num2str(lx)]);
                end
            elseif (B_type == 'y')
                for i = B_win
                    ori_signal                                      = squeeze(obj.img(index, i, :));
                    % wiener deconvolve
                    wienerdeconv_s                                  = fx_wiener_deconvolution_1d(ori_signal', kernel', q_factor);
                    img_WienerDeconv_temp(i - min(B_win) + 1, :)    = abs(hilbert(wienerdeconv_s));
                    % Wiener deconv. + ARregression
                    deconv                                          = class_deconv(ori_signal, kernel, obj.fs);
                    wienerdeconv_ARextrap_s                         = deconv.wienerdeconv_ARextrap(q_factor_AR, bw, k, 'bg', fft_padding);
                    img_WienerDeconv_AR_temp(i - min(B_win) + 1, :) = abs(hilbert(wienerdeconv_ARextrap_s(1:lz)));
%                     % sparse deconvolution
%                     % downsample the kernel
%                     kernel_sparse = kernel(1:DownRate:end);
%                     L = length(kernel_sparse);
%                     M = length(ori_signal);
%                     N = M + 1 - L;               % M : length of observed signal
%                     % Create convolution matrix H
%                     H = sparse(M, N);
%                     e = ones(N, 1);
%                     for ii = 0:L-1
%                         H = H + spdiags(kernel_sparse(ii + 1) * e, -ii, M, N);            % H : convolution matrix (sparse)
%                     end
%                     %                     issparse(H)                                       % confirm that H is a sparse matrix
%                     % L1 norm penalty The penalty function is phi(x) = lam abs(x)
%                     phi_L1 = @(x) lam * abs(x);
%                     wfun_L1 = @(x) abs(x)/lam;
%                     %                     dphi_L1 = @(x) lam * sign(x);
%                     [x1, ~] = fx_deconv_MM(ori_signal, phi_L1, wfun_L1, H, Nit);
%                     img_SparseDeconv_temp(i - min(B_win) + 1, :) = [zeros(round(L / 2 - 1), 1); x1; zeros(round(L / 2) - 1, 1)];
                    disp([num2str(i - min(B_win) + 1) '/' num2str(ly)]);
                end
            end
            obj.img_WienerDeconv_B    = img_WienerDeconv_temp;
            obj.img_WienerDeconv_AR_B = img_WienerDeconv_AR_temp;
            obj.img_SparseDeconv_B    = img_SparseDeconv_temp;
        end
        
        function obj = apply_deconvolutions(obj, q_factor, q_factor_AR, fft_padding, bw, bwr, k, ker_wid)
            % apply the deconvolutions on the original signals    
            % q_factor: the Q_factor for wiener deconvolution, 1e-2 for  default;
            % q_factor_AR: the Q_factor for wiener deconvolution combined with AR extrapolation;
            % fft_padding: the power of 2 for the length of the padded signal
            % bw: the bandwidht in unit dB
            % bwr: the the reference to define the bandwidht in unit dB
            % k: the order of AR regression
            % lam: the regularization parameter for sparse deconvolution
            % ker_wid: half of the width of the kernel, maximum centered
            % ****
            % move the maximum point to the center of the kernel
            kernel = obj.refer_Ascan;
            [~, I] = max(kernel);
            kernel = kernel(I - ker_wid: I + ker_wid);
            % remove the direct component
            kernel = kernel - mean(kernel);
            L      = length(kernel);
            n      = 2^nextpow2(L);
            Y      = fft(kernel, n);
            f      = obj.fs *(0:(n / 2)) / n;
            P      = abs(Y(1:n/2+1) / n);
            % make a series of f_windows
            f_windows                    = NaN(length(bwr), 2);
            for idx_bwr = 1:length(bwr)
                bwr_single               = bwr(idx_bwr);
                [f1, f2]                 = fx_calculate_bw(f, P, bwr_single);
                f_windows(idx_bwr, :)    = [f1, f2];
            end
            % create matrix
            [lx, ly, lz] = size(obj.img);
            % B scan image
            img_WienerDeconv_temp    = NaN(lx, ly, lz);
            img_WienerDeconv_AR_temp = NaN(lx, ly, lz);
            for i = 1:lx
                for j = 1:ly
                    ori_signal                        = squeeze(obj.img(i, j, :));
                    % remove the direct component
                    ori_signal                        = ori_signal - mean(ori_signal);
                    % wiener deconvolve
                    wienerdeconv_s                    = fx_wiener_deconvolution_1d(ori_signal', kernel', q_factor);
                    img_WienerDeconv_temp(i, j, :)    = abs(hilbert(wienerdeconv_s));
                    % Wiener deconv. + ARregression
                    deconv                            = class_deconv(ori_signal, kernel, obj.fs);
                    wienerdeconv_ARextrap_s           = deconv.wienerdeconv_ARextrap(q_factor_AR, f_windows, k, 'bg', fft_padding);
                    img_WienerDeconv_AR_temp(i, j, :) = abs(hilbert(wienerdeconv_ARextrap_s(1:lz)));
                end
                disp([num2str(i) '/' num2str(lx)]);
            end
            obj.img_WienerDeconv    = img_WienerDeconv_temp;
            obj.img_WienerDeconv_AR = img_WienerDeconv_AR_temp;
        end
        
        function obj = apply_deconvolutions_onlyWiener(obj, q_factor, ker_wid)
            % apply the deconvolutions on the original signals    
            % q_factor: the Q_factor for wiener deconvolution, 1e-2 for  default;
            % k: the order of AR regression
            % ker_wid: half of the width of the kernel, maximum centered
            % ****
            % move the maximum point to the center of the kernel
            kernel = obj.refer_Ascan;
            [~, I] = max(kernel);
            kernel = kernel(I - ker_wid: I + ker_wid);
            % remove the direct component
            kernel = kernel - mean(kernel);
            L      = length(kernel);
            n      = 2^nextpow2(L);
            Y      = fft(kernel, n);
            f      = obj.fs *(0:(n / 2)) / n;
            P      = abs(Y(1:n/2+1) / n);
            % create matrix
            [lx, ly, lz] = size(obj.img);
            % B scan image
            img_WienerDeconv_temp = NaN(lx, ly, lz);
            for i = 1:lx
                for j = 1:ly
                    ori_signal                     = squeeze(obj.img(i, j, :));
                    % remove the direct component
                    ori_signal                     = ori_signal - mean(ori_signal);
                    % wiener deconvolve
                    wienerdeconv_s                 = fx_wiener_deconvolution_1d(ori_signal', kernel', q_factor);
                    img_WienerDeconv_temp(i, j, :) = abs(hilbert(wienerdeconv_s));
                end
                clc;
                disp([num2str(i) '/' num2str(lx)]);
            end
            obj.img_WienerDeconv    = img_WienerDeconv_temp;
        end
        
        function demo_deconvolutions(obj, x, y, q_factor, q_factor_AR, fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit, fig_subname)
            % demonstrate origin signal and the deconvolved signal
            % q_factor: the Q_factor for wiener deconvolution, 1e-2 for  default;
            % q_factor_AR: the Q_factor for wiener deconvolution combined with AR extrapolation;
            % x, y: the x, y index to select the A scan
            % fft_padding: the power of 2 for the length of the padded signal
            % bw: the bandwidht in unit dB
            % bwr: the the reference to define the bandwidht in unit dB
            % k: the order of AR regression
            % lam: the regularization parameter for sparse deconvolution
            % ker_wid: half of the width of the kernel, maximum centered
            % DownRate: the parameter for downsampleing the kernel
            % Nit: the number of iteration for sparse deconvolution
            % fig_subname: the subname for the fig
            % ****
            % move the maximum point to the center of the kernel
            kernel = obj.refer_Ascan;
            [~, I] = max(kernel);
            kernel = kernel(I - ker_wid+1: I + ker_wid);
            L      = length(kernel);
            n      = 2^nextpow2(L);
            Y      = fft(kernel, n);
            f      = obj.fs *(0:(n / 2)) / n;
            P      = abs(Y(1:n/2+1) / n);
            % make a series of f_windows
            f_windows                    = NaN(length(bwr), 2);
            for idx_bwr = 1:length(bwr)
                bwr_single               = bwr(idx_bwr);
                [f1, f2]                 = fx_calculate_bw(f, P, bwr_single);
                f_windows(idx_bwr, :)    = [f1, f2];
            end
            figure('Name', ['Ref_signal_', fig_subname]);
            ca = subplot(2, 1, 1);
            plot(ca, (1:length(kernel)) / obj.fs * 1e6, kernel, 'linewidth', 2);
            ylabel(ca, '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            ca = subplot(2, 1, 2);   
            plot(ca, f / 1e6, P, 'linewidth', 2);
            ylabel(ca, '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel(ca, '\fontname {times new roman} Frequency (MHz)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            % plot origin signal
            ori_signal = squeeze(obj.img(x, y, :));
            % remove the direct component 
            ori_signal = ori_signal - mean(ori_signal);
            t_space    = (1:length(ori_signal)) / obj.fs * 1e6;
            cf = figure('Name', ['deconvolutions_', fig_subname]);
            set(cf, 'Position', [0, 0, 800, 1000], 'color', 'white');
            ca = subplot(3, 1, 1);
            plot (ca, t_space, ori_signal, 'linewidth', 2);
            hold on;
            plot (ca, t_space, abs(hilbert(ori_signal)), 'linewidth', 2);
            hold on;
            ylabel(ca, '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
%             xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            title('Original signal', 'fontsize', 16, 'Fontname', 'times new Roman');
            legend({'Recorded signal', 'Envelope of th recorded signal'});
            % wiener deconvolve
            deconvolved = fx_wiener_deconvolution_1d(ori_signal', kernel', q_factor);
            ca          = subplot(3, 1, 2);
            plot(ca, t_space, deconvolved, 'linewidth', 2);
            hold on;
            plot(ca, t_space, abs(hilbert(deconvolved)), 'linewidth', 2);
            legend({'Deconvolved signal', 'Envelope of the deconvolved signal'});
            ylabel(ca, '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
%             xlabel(ca, '\fontname {times new roman} time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            title('Wiener deconvolution', 'fontsize', 16, 'Fontname', 'times new Roman');
            % Wiener deconv. + ARregression
            deconv                      = class_deconv(ori_signal, kernel, obj.fs);
            wienerdeconv_ARextrap_s     = deconv.wienerdeconv_ARextrap(q_factor_AR, f_windows, k, 'bg', fft_padding);
            wienerdeconv_ARextrap_align = wienerdeconv_ARextrap_s(1:length(t_space));
            ca = subplot(3, 1, 3);
            plot(ca, t_space, wienerdeconv_ARextrap_align, 'linewidth', 2);
            hold on;
            plot(ca, t_space, abs(hilbert(wienerdeconv_ARextrap_align)), 'linewidth', 2);
            legend({'Deconvolved signal', 'Envelope of the deconvolved signal'});
            ylabel(ca, '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
            xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
%             ylim([min(wienerdeconv_ARextrap_s) / 5 max(wienerdeconv_ARextrap_s) / 5]);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            title('Wiener deconvolution combined with AR extrapolation', 'fontsize', 16, 'Fontname', 'times new Roman');
%             % sparse deconvolution
%             % downsample the kernel
%             kernel_sparse = kernel(1:DownRate:end);
%             L = length(kernel_sparse);
%             M = length(ori_signal);
%             N = M + 1 - L;               % M : length of observed signal
%             % Create convolution matrix H
%             H = sparse(M, N);
%             e = ones(N, 1);
%             for ii = 0:L-1
%                 H = H + spdiags(kernel_sparse(ii + 1) * e, -ii, M, N);            % H : convolution matrix (sparse)
%             end
%             issparse(H)                                       % confirm that H is a sparse matrix
%             % L1 norm penalty The penalty function is phi(x) = lam abs(x)
%             phi_L1 = @(x) lam * abs(x);
%             wfun_L1 = @(x) abs(x)/lam;
% %             dphi_L1 = @(x) lam * sign(x);
%             [x1, cost1] = fx_deconv_MM(ori_signal, phi_L1, wfun_L1, H, Nit);
%             recons_signal = [zeros(round(L / 2), 1); x1; zeros(round(L / 2), 1)];
%             ca = subplot(4, 1, 4);
%             plot(ca, t_space, recons_signal(1:length(t_space)), 'linewidth', 2);
%             ylabel(ca, '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
%             xlabel(ca, '\fontname {times new roman} Time(\mus)', 'fontsize', 16);
% %             ylim([min(wienerdeconv_ARextrap_s) / 5 max(recons_signal) / 5]);
%             set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
%             set(ca, 'linewidth', 2);
%             title('Sparse deconvolution', 'fontsize', 16, 'Fontname', 'times new Roman');
%             cf = figure();
%             plot(cost1);
        end
    
        function demo_deconvolutions_Bscans(obj, B_type, index, Bwin, fig_subname)
            % demonstrate the B_scan of the deconvolutions resutls
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % Bwin: the indexes of the first and the last points in Bscan 
            % image
            % fig_subname: the subname for the fig
            img_WienerDeconv_temp    = obj.img_WienerDeconv_B(Bwin - min(Bwin) + 1, :);
            img_WienerDeconv_AR_temp = obj.img_WienerDeconv_AR_B(Bwin - min(Bwin) + 1, :);
            img_SparseDeconv_temp    = obj.img_SparseDeconv_B(Bwin - min(Bwin) + 1, :); 
            x                        = (1: size(img_WienerDeconv_temp, 1) ) / obj.fx * 1e3;
            z                        = (1: size(img_WienerDeconv_temp, 2) ) / obj.fs * 1e6;
            cf = figure('Name', ['Bscan_WienerDeconv' '_', B_type, '_', num2str(index), fig_subname]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            imagesc(img_WienerDeconv_temp');
            % figure setup
            colormap(gray);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            cl   = caxis;
            caxis(cl / 5);
%             caxis([min(img_WienerDeconv_temp, [], 'all') / 5 max(img_WienerDeconv_temp, [], 'all') / 10 ]);
            %
            cf = figure('Name', ['Bscan_WienerDeconv_AR' '_', B_type, '_', num2str(index), fig_subname]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            imagesc(x, z, img_WienerDeconv_AR_temp');
            % figure setup
            colormap(jet);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            caxis('auto');
            cl   = caxis;
            caxis(cl / 5);
            %
            colormap(gray);
            h  = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Magnitude (arb.)\newline ');
%             cf = figure('Name', ['Bscan_SparseDeconv' '_', B_type, '_', num2str(index), fig_subname]);
%             set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
%             imagesc(x, z, img_SparseDeconv_temp');
%             % figure setup
%             colormap(gray);
%             h  = colorbar;
%             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
%             set(gca, 'Fontname', 'times new Roman');
%             xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
%             ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
%             set(gca, 'fontsize', 16);
%             set(gca, 'linewidth', 2);
%             caxis('auto');
%             cl =  caxis;
%             caxis(cl / 10);
        end
        
        function obj = track_interply_inam(obj, PropertyName, nol)
            % track the interply by instantaneous amp..
            % return a masked dataset
            % PropertyName: 'img_WienerDeconv' or 'img_WienerDeconv_AR' ...
            % nol: number of layers to track
            % ************
            % The sequence is changed unexpectedly. Transverse the matrix
            img_ori            = obj.(PropertyName);
            mask_plytrack_temp = zeros(size(img_ori));
            row_temp           = [];
            col_temp           = [];
            dep_temp           = [];
            layer_counts_temp  = [];
            threshold          = 0.0; % useless now
            errors_pos_repeat  = NaN(size(img_ori, 1),size(img_ori, 2), nol);
            lct_3d             = NaN(size(img_ori, 1), size(img_ori, 2), nol);
            thickness          = NaN(size(mask_plytrack_temp,1), ...
                size(mask_plytrack_temp,2), nol-1);
            for idx = 1:size(img_ori, 1)
                for idy = 1:size(img_ori, 2)
                    if obj.rear_I(idx, idy) <=0 || obj.rear_I(idx, idy) > size(img_ori, 3) ...
                            || obj.front_I(idx, idy) <=0 || obj.front_I(idx, idy) > obj.rear_I(idx, idy)
                        continue;
                    end
                    A_ori                  = abs(squeeze(img_ori(idx, idy, :)));
                    % track the phase by wrapped inst. phase
                    %                     Idx_resin                                  = linspace(obj.front_I(idx, idy), obj.rear_I(idx, idy), 25);
                    %                     [lct, ~, ~, ~, ~]                          = fx_track_byphase(Idx_resin, A_ori, threshold);
                    % track the -pi/2 by distance
                    %                     [lct, ~, ~] = fx_track_bycluster(hilbert(A_ori), threshold, 23);
                    Idx_resin              = linspace(obj.front_I(idx, idy), ...
                        obj.rear_I(idx, idy), 25);  % !!!!! 25 is the number of layers incl. front and back surfaces
                    % track the -pi/2 by phase
                    power_noise            = 1e-10;
                    [lct, ~, er_pos, ~]    = fx_track_bypeaks(Idx_resin, real(A_ori), threshold, power_noise);
                    % save to 3d vector
                    lct_3d(idx, idy, :)    = lct;
                    %                     mask_plytrack_temp(idx, idy, lct(2:end-1)) = 1;
                    thickness(idx, idy, :) = diff(lct);
                    % derive the rear surface by phase
                    %                     I_rear_surface       = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                    %                         phase_rear_surface) + round(obj.front_I(idx, idy));
                    %                     obj.rear_I(idx, idy)   = lct(end);
                    obj.front_I(idx, idy)  = lct(1);
                    % save the 3D scatter points
                    num_interplies         = length(lct)-2; % incl. the front and rear surfaces
                    layer_counts_temp      = cat(2, layer_counts_temp, (1:num_interplies));
                    row_temp               = cat(2, row_temp, idx * ones(1, num_interplies));
                    col_temp               = cat(2, col_temp, idy * ones(1, num_interplies));
                    dep_temp               = cat(2, dep_temp, round(lct(2:end-1)));
                    errors_pos_repeat(idx, idy, :) = er_pos;
                end
                disp([num2str(idx), '/' , num2str(size(img_ori, 1))]);
            end
            lct_3d_interp = round(lct_3d);
            % interpolation
            %             for i = 1:25
            %                 lct_3d_interp(:, :, i) = round(fx_inpaint_nans(lct_3d(:, :, i), 1));
            %             end
            % fill mask vector
            for idx = 1:size(img_ori, 1)
                for idy = 1:size(img_ori, 2)
                    % check nan
                    for idz = 2:size(lct_3d_interp, 3)-1
                        if ~isnan(lct_3d_interp(idx, idy, idz))
                            mask_plytrack_temp(idx, idy, lct_3d_interp(idx, idy, idz)) = 1;
                        end
                    end
                end
            end
            obj.mask_plytrack = mask_plytrack_temp;
            obj.row           = row_temp;
            obj.col           = col_temp;
            obj.dep           = dep_temp;
            obj.layer_counts  = layer_counts_temp;
            obj.errors_pos    = errors_pos_repeat;
            obj.est_thickness = thickness;
            % error bar
            error_mean = squeeze(mean(abs(errors_pos_repeat),[1 2], 'omitnan'));
            error_std  = squeeze(std(errors_pos_repeat, 0, [1 2], 'omitnan'));
            cf = figure('Name', ['error_bar' '_' PropertyName]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(1, 1, 1);
            errorbar((1:nol)/obj.fs*1e6, error_mean*100, error_std*100,'-s',...
                'MarkerSize', 10, 'Linewidth', 3,'MarkerEdgeColor','red','MarkerFaceColor','red')
            xlabel('\fontname{times new roman} Time (µm)', 'fontsize', 16);
            ylabel('\fontname{times new roman} Errors (mean and std)', 'fontsize', 16);
            ytickformat('percentage');
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            grid on;
        end
        
        function obj = track_interply_inam_BFsurfacesbased(obj, PropertyName, nol)
            % track the interply by instantaneous amp..
            % return a masked dataset
            % PropertyName: 'img_WienerDeconv' or 'img_WienerDeconv_AR' ...
            % nol: number of layers to track
            % ************
            % The sequence is changed unexpectedly. Transverse the matrix
            img_ori            = obj.(PropertyName);
            mask_plytrack_temp = zeros(size(img_ori));
            row_temp           = [];
            col_temp           = [];
            dep_temp           = [];
            layer_counts_temp  = [];
            threshold          = 0.0; % useless now
            errors_pos_repeat  = NaN(size(img_ori, 1),size(img_ori, 2), nol);
            lct_3d             = NaN(size(img_ori, 1), size(img_ori, 2), nol);
            thickness          = NaN(size(mask_plytrack_temp,1), ...
                size(mask_plytrack_temp,2), nol-1);
            for idx = 1:size(img_ori, 1)
                for idy = 1:size(img_ori, 2)
                    if obj.rear_I(idx, idy) <=0 || obj.rear_I(idx, idy) > size(img_ori, 3) ...
                            || obj.front_I(idx, idy) <=0 || obj.front_I(idx, idy) > obj.rear_I(idx, idy)
                        continue;
                    end
                    A_ori                  = squeeze(img_ori(idx, idy, :));
                    % track the phase by wrapped inst. phase
                    %                     Idx_resin                                  = linspace(obj.front_I(idx, idy), obj.rear_I(idx, idy), 25);
                    %                     [lct, ~, ~, ~, ~]                          = fx_track_byphase(Idx_resin, A_ori, threshold);
                    % track the -pi/2 by distance
                    %                     [lct, ~, ~] = fx_track_bycluster(hilbert(A_ori), threshold, 23);
                    Idx_resin              = linspace(obj.front_I(idx, idy), ...
                        obj.rear_I(idx, idy), 25);  % !!!!! 25 is the number of layers incl. front and back surfaces
                    % track the -pi/2 by phase
                    power_noise            = 1e-10;
                    [lct, ~, er_pos, ~]    = fx_track_bypeaks_BFsurfacesbased(Idx_resin, real(A_ori), threshold, power_noise);
                    % save to 3d vector
                    lct_3d(idx, idy, :)    = lct;
                    %                     mask_plytrack_temp(idx, idy, lct(2:end-1)) = 1;
                    thickness(idx, idy, :) = diff(lct);
                    % derive the rear surface by phase
                    %                     I_rear_surface       = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                    %                         phase_rear_surface) + round(obj.front_I(idx, idy));
                    %                     obj.rear_I(idx, idy)   = lct(end);
                    obj.front_I(idx, idy)  = lct(1);
                    % save the 3D scatter points
                    num_interplies         = length(lct)-2; % incl. the front and rear surfaces
                    layer_counts_temp      = cat(2, layer_counts_temp, (1:num_interplies));
                    row_temp               = cat(2, row_temp, idx * ones(1, num_interplies));
                    col_temp               = cat(2, col_temp, idy * ones(1, num_interplies));
                    dep_temp               = cat(2, dep_temp, round(lct(2:end-1)));
                    errors_pos_repeat(idx, idy, :) = er_pos;
                end
                disp([num2str(idx), '/' , num2str(size(img_ori, 1))]);
            end
            lct_3d_interp = round(lct_3d);
            % interpolation
            %             for i = 1:25
            %                 lct_3d_interp(:, :, i) = round(fx_inpaint_nans(lct_3d(:, :, i), 1));
            %             end
            % fill mask vector
            for idx = 1:size(img_ori, 1)
                for idy = 1:size(img_ori, 2)
                    % check nan
                    for idz = 2:size(lct_3d_interp, 3)-1
                        if ~isnan(lct_3d_interp(idx, idy, idz))
                            mask_plytrack_temp(idx, idy, lct_3d_interp(idx, idy, idz)) = 1;
                        end
                    end
                end
            end
            obj.mask_plytrack = mask_plytrack_temp;
            obj.row           = row_temp;
            obj.col           = col_temp;
            obj.dep           = dep_temp;
            obj.layer_counts  = layer_counts_temp;
            obj.errors_pos    = errors_pos_repeat;
            obj.est_thickness = thickness;
            % error bar
            error_mean = squeeze(mean(abs(errors_pos_repeat),[1 2], 'omitnan'));
            error_std  = squeeze(std(errors_pos_repeat, 0, [1 2], 'omitnan'));
            cf = figure('Name', ['error_bar' '_' PropertyName]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(1, 1, 1);
            errorbar((1:nol)/obj.fs*1e6, error_mean*100, error_std*100,'-s',...
                'MarkerSize', 10, 'Linewidth', 3,'MarkerEdgeColor','red','MarkerFaceColor','red')
            xlabel('\fontname{times new roman} Time (µm)', 'fontsize', 16);
            ylabel('\fontname{times new roman} Errors (mean and std)', 'fontsize', 16);
            ytickformat('percentage');
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ca, 'linewidth', 2);
            grid on;
        end
        
        % ************************* In-plane direction RT *********
        function [obj, C_scan_inam_temp, C_scan_index_temp] = define_parallel_inamCscan(obj, ratio, PropertyName)
            % demonstrate the profile of one inter-ply in 3d
            % save the C_scan_inam
            % ratio: the ratio determing the index by the front and back surfaces , unit: arb.
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % ***
            if strcmp(PropertyName, 'img')
                inam = obj.(PropertyName);
            else
                inam = abs(obj.(PropertyName));
            end
            C_scan_inam_temp  = NaN([size(inam, 1), size(inam, 2)]);
            C_scan_index_temp = NaN([size(inam, 1), size(inam, 2)]);
            for i = 1: size(inam, 1)
                for j = 1:size(inam, 2)
                    % take both front and back surfaces into consideration
                    front_index             = obj.front_I(i, j);
                    back_index              = obj.rear_I(i, j);
                    if length(ratio)>1 % take average of the amplitude
                        C_scan_inam_temp(i, j) = 0;
                        for r_index = 1:length(ratio) % loop the depth 
                            index                   = (1 - ratio(r_index))*front_index + ratio(r_index)*back_index;
%                             C_scan_index_temp(i, j) = index;
                            if index < back_index && ~isnan(index)
                                C_scan_inam_temp(i, j) = C_scan_inam_temp(i, j) + inam(i, j, round(index));
                            end
                        end
                    else % no avearage
                        index                   = (1 - ratio)*front_index + ratio*back_index;
                        C_scan_index_temp(i, j) = index;
                        if index < back_index && ~isnan(index)
                            C_scan_inam_temp(i, j) = inam(i, j, round(index));
                        end
                    end
                end
            end
            % filter the profile and asign the inam
            C_scan_inam_temp = medfilt2(C_scan_inam_temp);
            C_scan_inam_temp = fx_inpaint_nans(C_scan_inam_temp, 5);
            obj.C_scan_inam  = C_scan_inam_temp;
%             %
%             figure('Name', ['Cscan_RT_parallel_' , PropertyName, '_', num2str(ratio)]);
%             ax = subplot(1, 1, 1);
%             X = (1: size(obj.C_scan_inam, 2)) / obj.fx * 1e3;
%             Y = (1: size(obj.C_scan_inam, 1)) / obj.fy * 1e3;
%             imagesc(ax, X, Y, obj.C_scan_inam);
%             hold on;
%             h = colorbar;
%             colormap jet;
%             %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
%             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} inst. amp. (arb,)');
%             set(ax, 'fontsize', 16);
%             set(ax, 'linewidth', 1.5);
%             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
%             ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
%             %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
%             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
%             set(gca, 'linewidth', 2);
        end
        
        function [obj, C_scan_inam_ply, C_scan_index] = define_plywise_inamCscan(obj, layer, ratios, PropertyName)
            % demonstrate the profile of one inter-ply in 3d
            % save the C_scan_inam
            % layer: number of the layer
            % ratios: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % ***
            inam             = abs(obj.(PropertyName));
            inph_ex          = obj.mask_plytrack;
            C_scan_index     = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_inam_temp = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_inam_ply  = zeros([size(inph_ex, 1), size(inph_ex, 2)]);
            for r = 1: length(ratios)
                ratio = ratios(r);
                for i = 1: size(inph_ex, 1)
                    for j = 1:size(inph_ex, 2)
                        if round(obj.front_I(i, j))
                            inph_ex(i, j, round(obj.front_I(i, j))) = 1;
                        end
                        % no back_surface or delamination
                        if ~isnan(obj.rear_I(i, j))
                            inph_ex(i, j, round(obj.rear_I(i, j))) = 1;
                        end
                        k = find(inph_ex(i, j, :), layer + 1);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                        % extract two inter-plies embedding the ply
                        if length(k) > layer
                            k1                 = k(layer);
                            k2                 = k(layer + 1);
                            % use the ratio to calculate the position in the ply
                            z_index            = round(k1 + (k2 - k1) * ratio);
                            C_scan_index(i, j) = z_index;
                        else % set the last inter-ply as NaN otherwise
                            C_scan_index(i, j) = NaN;
                        end
                    end
                    
                end
                % Perform a morphological close operation on the image.
                NaN_part = ~isnan(C_scan_index);
                se       = strel('disk', 5);
                closeNaN = ~imclose(~NaN_part,se);    
                % fillnan
                C_scan_index = fx_inpaint_nans(C_scan_index, 5);
%                 C_scan_index(~closeNaN) = NaN;
                % filter the profile and asign the inam
                C_scan_index = medfilt2(C_scan_index, [5 5], 'symmetric');
                C_scan_index(~closeNaN) = NaN;
%                 % limit the variation
%                 C_scan_index_ave = mean(C_scan_index, 'all', 'omitnan');
%                 C_scan_index(C_scan_index<C_scan_index_ave-19) = C_scan_index_ave;
%                 C_scan_index(C_scan_index>C_scan_index_ave+19) = C_scan_index_ave;
                C_scan_index     = round(C_scan_index);
                for i = 1: size(inph_ex, 1)
                    for j = 1:size(inph_ex, 2)
                        if obj.rear_I(i, j) > C_scan_index(i, j) && inam(i, j, C_scan_index(i, j)) % check the rear surface and assign NaN
                            C_scan_inam_temp(i, j) = inam(i, j, C_scan_index(i, j));
                        else
                            C_scan_inam_temp(i, j) = NaN;
                        end
                    end
                end
                C_scan_inam_ply = C_scan_inam_ply + C_scan_inam_temp;
            end
            C_scan_inam_ply = C_scan_inam_ply / length(ratios);
%             % Perform a morphological close operation on the image.
%             NaN_part        = ~isnan(C_scan_inam_ply);
%             se              = strel('disk', 5);
%             closeNaN        = imclose(NaN_part,se);
%             % fillnan
%             C_scan_inam_ply            = fx_inpaint_nans(C_scan_inam_ply, 5);
%             C_scan_inam_ply(~closeNaN) = NaN;
%             C_scan_inam_ply            = medfilt2(C_scan_inam_ply, [3 3], 'symmetric');
          % commented out for waviness sample evaluation
%             check the variation to find the large echo from the delamination
            inam_mean = mean(C_scan_inam_ply, [1 2], 'omitnan');
            inam_std  = std(C_scan_inam_ply, [], 'all', 'omitnan');
            C_scan_inam_ply(C_scan_inam_ply>inam_mean+6*inam_std) = NaN;
            obj.C_scan_inam_plywise = C_scan_inam_ply;
%             % imshow
%             figure; clf; imagesc(C_scan_inam_ply); colormap gray; axis('image');
%             title('ply-wise C-scan Image');
        end
        
        function obj = show_one_Radontransform(obj, layer, ratio, PropertyName, center, radii)
            % demonstrate the profile of one inter-ply in 3d
            % save the C_scan_inam
            % layer: number of the layer
            % ratio: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % ***
            % the ratio should be 0 ~ 1
            if ratio < 0 || ratio > 1
                msg = "The 'raio' should be between 0 and 1";
                error(msg);
            end
            inam            = abs(obj.(PropertyName));
            inph_ex         = obj.mask_plytrack;
            profile_layer   = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_index    = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            obj.C_scan_inam = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    if round(obj.front_I(i, j))
                        inph_ex(i, j, round(obj.front_I(i, j))) = 1;
                    end
                    if round(obj.rear_I(i, j))
                        inph_ex(i, j, round(obj.rear_I(i, j))) = 1;
                    end
                    k = find(inph_ex(i, j, :), layer + 1);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                    % extract two inter-plies embedding the ply
                    if length(k) > layer
                        k1 = k(layer);
                        k2 = k(layer + 1);
                        % use the ratio to calculate the position in the ply
                        z_index = round(k1 + (k2 - k1) * ratio);
                        C_scan_index(i, j)    = z_index;
                        obj.C_scan_inam(i, j) = inam(i, j, z_index);
                    else % set the last inter-ply as NaN otherwise
                        C_scan_index(i, j)    = NaN;
                    end
                end
            end
            % filter the profile and asign the inam
            C_scan_index = medfilt2(C_scan_index);
            C_scan_index = fx_inpaint_nans(C_scan_index, 5);
            % average the value exceed the range
            C_scan_index_ave = round(mean(C_scan_index, [1 2], 'omitnan'));
            C_scan_index(C_scan_index<C_scan_index_ave-19) = C_scan_index_ave;
            C_scan_index(C_scan_index>C_scan_index_ave+19) = C_scan_index_ave;
            profile_layer = C_scan_index / obj.fs * 1e6;
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    obj.C_scan_inam(i, j) = inam(i, j, C_scan_index(i, j));
                end
            end
            %             % fillna
%             obj.C_scan_inam(isnan(obj.C_scan_inam)) = mean(obj.C_scan_inam, 'all', 'omitnan');
            % Radon transform
            % creat a round mask parameters
            theta                            = 0:5:180;
            [R, anguler_1D, max_angle_I, xp] = fx_Radonto1Dangular(obj.C_scan_inam, center, radii, theta);
            % plot the R-theta image
            figure('Name', ['Radon_transform_R_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, theta, xp, R);
            hold on;
            h = colorbar;
            % plot the 1D anguler (sum of gradient by central-different kernel)
            figure('Name', ['Radon_transform_1Danguler_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            plot(theta, anguler_1D, 'linewidth', 2);
            % creat a round mask
            mask = fx_createCirclesMask(obj.C_scan_inam, center, radii);
            x1 = center(2) - radii;
            x2 = center(2) + radii;
            y1 = center(1) - radii;
            y2 = center(1) + radii;
            % use the mask to select a round region
            C_scan_inam_mask = obj.C_scan_inam .* mask;
            C_scan_inam_mask = C_scan_inam_mask(x1: x2,y1: y2);
            figure('Name', ['Cscan_RT_fitplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio)]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_mask, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan_inam_mask, 2)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_mask);
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} inst. amp. (arb,)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % plot the profile of the layer
            profile_layer = profile_layer .* mask;
            profile_layer = profile_layer(x1: x2, y1: y2);
            figure('Name', ['profile_RT_fitplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio)]);
            ax = subplot(1, 1, 1);
            X = (1: size(profile_layer, 1)) / obj.fx * 1e3;
            Y = (1: size(profile_layer, 2)) / obj.fy * 1e3;
            imagesc(ax, X,Y, profile_layer)
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function [max_angle_I, obj]= show_one_Radontransform_vitualprofile(obj, layer, ratio, PropertyName, center, radii, max_layer)
            % demonstrate the profile of one inter-ply in 3d
            % save the C_scan_inam
            % layer: number of the layer, 'layer'th layer
            % ratio: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % max_layer: the total plies of the structure
            % surfaces.
            % ***
            % the ratio should be 0 ~ 1
            if ratio < 0 || ratio > 1
                msg = "The 'raio' should be between 0 and 1";
                error(msg);
            end
            inam = abs(obj.(PropertyName));
            % 'layer'th interply
            oneply_row = obj.row(obj.layer_counts_new==layer);
            oneply_col = obj.col(obj.layer_counts_new==layer);  % / obj.fy * 1e3 ;
            oneply_dep = obj.dep(obj.layer_counts_new==layer);  % / obj.fy * 1e3 ;
            % 'layer - 1'th interply
            oneply_row_pre = obj.row(obj.layer_counts_new==layer - 1);
            oneply_col_pre = obj.col(obj.layer_counts_new==layer - 1);  % / obj.fy * 1e3 ;
            oneply_dep_pre = obj.dep(obj.layer_counts_new==layer - 1);  % / obj.fy * 1e3 ;
            % construct the profile of 'layer'th and 'layer - 1'th
            if layer == 1 % front surface
                profile_layer_pre = obj.front_I;
            else
                profile_layer_pre = NaN(max(oneply_row_pre), max(oneply_col_pre));
                profile_layer_z = NaN(max(oneply_row_pre), max(oneply_col_pre));
                for i = 1: length(oneply_row_pre)
                    profile_layer_pre(oneply_row_pre(i), oneply_col_pre(i)) =  oneply_dep_pre(i);
                end
            end
            if layer == max_layer % rear surface
                profile_layer = obj.rear_I;
            else
                profile_layer = NaN(max(oneply_row), max(oneply_col));
                profile_layer_z = NaN(max(oneply_row), max(oneply_col));
                for i = 1: length(oneply_row)
                    profile_layer(oneply_row(i), oneply_col(i)) =  oneply_dep(i);
                end
            end
            % construct the target profile and the inam
            obj.C_scan_inam = NaN(size(profile_layer));
            profile_layer_ave = mean(profile_layer, 'all', 'omitnan');
            profile_layer_pre_ave = mean(profile_layer_pre, 'all', 'omitnan');
            %
            for i = 1:size(profile_layer, 1)
                for j = 1:size(profile_layer, 2)
                    if isnan(profile_layer_pre(i, j ))
                        k1 = profile_layer_pre_ave;
                    else
                        k1 = profile_layer_pre(i, j );
                    end
                    if isnan(profile_layer(i, j ))
                        k2 = profile_layer_ave;
                    else
                        k2 = profile_layer(i, j );
                    end
                    if k2 > k1 && obj.rear_I(i, j) >= k2
                        % layer is deeper than previous layer & rear_surface is deeper than all
                        z_index = round(k1 + (k2 - k1) * ratio);
                        obj.C_scan_inam(i, j) = inam(i, j, z_index);
                        profile_layer_z(i, j) = z_index/ obj.fs * 1e6;
                    end
                end
            end
            % fillna
            obj.C_scan_inam(isnan(obj.C_scan_inam)) = mean(obj.C_scan_inam, 'all', 'omitnan');
            % Radon transform
            % creat a round mask parameters
            [R, anguler_1D, max_angle_I, xp, theta] = fx_Radonto1Dangular(obj.C_scan_inam, center, radii);
            % plot the R-theta image
            figure('Name', ['Radon_transform_R_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, theta, xp, R);
            hold on;
            h = colorbar;
            % plot the 1D anguler (sum of gradient by central-different kernel)
            figure('Name', ['Radon_transform_1Danguler_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            plot(theta, anguler_1D, 'linewidth', 2);
            % creat a round mask
            mask = fx_createCirclesMask(obj.C_scan_inam, center, radii);
            x1 = center(2) - radii;
            x2 = center(2) + radii;
            y1 = center(1) - radii;
            y2 = center(1) + radii;
            % use the mask to select a round region
            C_scan_inam_mask = obj.C_scan_inam .* mask;
            C_scan_inam_mask = C_scan_inam_mask(x1: x2,y1: y2);
            figure('Name', ['Cscan_RT_fitplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio)]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_mask, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan_inam_mask, 2)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_mask);
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} inst. amp. (arb,)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % plot the profile of the layer
            profile_layer_z = profile_layer_z .* mask;
            profile_layer_z = profile_layer_z(x1: x2, y1: y2);
            figure('Name', ['profile_RT_fitplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio)]);
            ax = subplot(1, 1, 1);
            X = (1: size(profile_layer_z, 1)) / obj.fx * 1e3;
            Y = (1: size(profile_layer_z, 2)) / obj.fy * 1e3;
            imagesc(ax, X,Y, profile_layer_z)
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function [Cscan_inam, profile_layer, profile_layer_pre] = define_Cscan_depthprofile_knn(obj, layer, ratio, PropertyName, max_layer)
            % define the profile of one inter-ply in 3d
            % save the C_scan_inam
            % layer: number of the layer, 'layer'th layer
            % ratio: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % max_layer: the total plies of the structure
            % surfaces.
            % ***
            % the ratio should be 0 ~ 1
            if ratio < 0 || ratio > 1
                msg = "The 'raio' should be between 0 and 1";
                error(msg);
            end
            inam = abs(obj.(PropertyName));
            if isempty(obj.layer_counts_new) % no knn and no layer_counts_new
                layer_counts_temp = obj.layer_counts;
            else
                layer_counts_temp = obj.layer_counts_new;
            end
            % 'layer'th interply
            oneply_row = obj.row(layer_counts_temp==layer);
            oneply_col = obj.col(layer_counts_temp==layer);  % / obj.fy * 1e3 ;
            oneply_dep = obj.dep(layer_counts_temp==layer);  % / obj.fy * 1e3 ;
            % 'layer - 1'th interply
            oneply_row_pre = obj.row(layer_counts_temp==layer - 1);
            oneply_col_pre = obj.col(layer_counts_temp==layer - 1);  % / obj.fy * 1e3 ;
            oneply_dep_pre = obj.dep(layer_counts_temp==layer - 1);  % / obj.fy * 1e3 ;
            % construct the profile of 'layer'th and 'layer - 1'th
            if layer == 1 % front surface
                profile_layer_pre = obj.front_I;
            else
                profile_layer_pre = NaN(max(oneply_row_pre), max(oneply_col_pre));
                profile_layer_z   = NaN(max(oneply_row_pre), max(oneply_col_pre));
                for i = 1: length(oneply_row_pre)
                    profile_layer_pre(oneply_row_pre(i), oneply_col_pre(i)) =  oneply_dep_pre(i);
                end
            end
            if layer == max_layer % rear surface
                profile_layer = obj.rear_I;
            else
                profile_layer   = NaN(max(oneply_row), max(oneply_col));
                profile_layer_z = NaN(max(oneply_row), max(oneply_col));
                for i = 1: length(oneply_row)
                    profile_layer(oneply_row(i), oneply_col(i)) =  oneply_dep(i);
                end
            end
            % construct the target profile and the inam
            Cscan_inam            = NaN(size(profile_layer));
            profile_layer_ave     = mean(profile_layer, 'all', 'omitnan');
            profile_layer_pre_ave = mean(profile_layer_pre, 'all', 'omitnan');
            %
            for i = 1:size(profile_layer, 1)
                for j = 1:size(profile_layer, 2)
                    if i <= size(profile_layer_pre, 1) &&... % avoid exceeding the bound
                        j <= size(profile_layer_pre, 2) &&...
                        isnan(profile_layer_pre(i, j))
                        k1 = profile_layer_pre_ave;
                    elseif i <= size(profile_layer_pre, 1) &&... % avoid exceeding the bound
                        j <= size(profile_layer_pre, 2)
                        k1 = profile_layer_pre(i, j );
                    else
                        k1 = profile_layer_pre_ave;
                    end
                    if isnan(profile_layer(i, j ))
                        k2 = profile_layer_ave;
                    else
                        k2 = profile_layer(i, j );
                    end
                    if k2 > k1 && obj.rear_I(i, j) >= k2 % k1 exists and smaller than k2 
                        % layer is deeper than previous layer & rear_surface is deeper than all
                        z_index               = round(k1 + (k2 - k1) * ratio);
                        Cscan_inam(i, j)      = inam(i, j, z_index);
                        profile_layer_z(i, j) = z_index;
                    end
                end
            end
            % fillna
            Cscan_inam(isnan(Cscan_inam)) = mean(Cscan_inam, 'all', 'omitnan');
        end
        
        function obj = compute_orientation_by_RT_correct(obj, radius, theta, imagename)
            % extract the inplane orientation by Radon transform
            % radius: the radius array for RT
            % theta: angle space
            % imagename: the name of the 2d image
            % ***
            % fillna
            C_scan_inam_RT           = obj.(imagename);
            NaN_part                 = isnan(C_scan_inam_RT);
%             % fill nan
%             C_scan_inam_RT           = fx_inpaint_nans(C_scan_inam_RT, 5);
            %             C_scan_inam_RT(isnan(C_scan_inam_RT)) = mean(C_scan_inam_RT, 'all', 'omitnan');
            [lxc, lyc]               = size(C_scan_inam_RT);
            ID                       = NaN(lxc, lyc, length(radius), length(theta));
            %             ID_gyy = NaN(lxc, lyc, length(radiis), length(theta));
            % pad the C_scan_inam matrix
            r_max                    = max(radius);
            C_scan_inam_RT_pad       = C_scan_inam_RT;
            %             C_scan_inam_RT_pad       = padarray(C_scan_inam_RT, [r_max r_max], mean(C_scan_inam_RT, 'all', 'omitnan'));
            % check the nan and adjust the radii
            [r_nan, l_nan] = find(NaN_part);
            % timer
            tic;
            for i = 1+r_max: lxc - r_max
                for j = 1+r_max: lyc - r_max
                    if (NaN_part(i, j))
                        continue;
                    end
                    for r_idx = 1: length(radius)
                        r          = radius(r_idx);
%                         check the nan and adjust the radii
                        distances = sqrt((r_nan-i).^2 + (l_nan-j).^2);
                        if r > min(distances)
                            r = floor(min(distances));
                        end
                        if r < 4
                            break;
                        end
                        center                = [j, i];
                        % creat a round mask parameters
                        [~, anguler_1D, ~, ~] = fx_Radonto1Dangular_correct(C_scan_inam_RT_pad, center, r, theta);
%                         [~, anguler_1D, ~, ~] = fx_Radonto1Dangular(C_scan_inam_RT_pad, center, r, theta);
                        ID(i, j, r_idx, :)    = anguler_1D;
                        %                         ID_gyy(i, j, r_idx, :) = anguler_1D_gyy;
                    end
                end
                clc;
                fprintf('RT progress: %0.2f%%\n',100*(i-r_max)/(lxc-r_max));
            end
            timeElapsed = toc;
            disp(['Compute_MultiResolution_RT ', num2str(timeElapsed)]);
            obj.RT_ID = ID; 
        end
        
        function show_ID_RT(obj, x, y, radius, theta)
            % demonstrate the information diagram of one point in the logGabor filtered images
            % use the logGabor filtered images calcualated before
            % x, y: index of the point in the image unit : pixel
            % radius, theta: the same ratio for the multi-resolution RT
            % ***
            r_max = max(radius);
            ID = squeeze(obj.RT_ID(x, y, :, :));
            %             ID_gyy = squeeze(obj.RT_ID_gyy(x, y, :, :));
            angular_1D = sum(ID, 1);
            [~, indAng] = max(angular_1D);
            disp(theta(indAng));
            % fillna
            C_scan_inam_LG = obj.C_scan_inam;
            C_scan_inam_LG(isnan(C_scan_inam_LG)) = mean(C_scan_inam_LG, 'all', 'omitnan');
            %
            figure('Name', ['Cscan_fitplytrack_forlogGabor', '_']);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_LG, 2)) / obj.fx * 1e3; % default X = (1: size(C_scan_inam_LG, 2))  here
            Y = (1: size(C_scan_inam_LG, 1)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_LG);
            hold on;
            scatter(ax, y / obj.fx * 1e3, x / obj.fy * 1e3, 30, 'r', 'filled'); % reversed x, y here
            colormap jet;
            h = colorbar;
            % caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            % zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % demonstrate the ID and 1D angular
            if length(radius) > 1
                cf = figure();
                ax = subplot(1, 1, 1);
                surf(theta, radius, ID);
                h = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
                set(ax, 'fontsize', 16);
                set(ax, 'linewidth', 1.5);
                xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
                ylabel('\fontname {times new roman} Wavelength (pixel)', 'fontsize', 16);
                %                 cf = figure();
                %                 ax = subplot(1, 1, 1);
                %                 surf(theta, radius, ID_gyy);
                %                 h = colorbar;
                %                 set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
                %                 set(ax, 'fontsize', 16);
                %                 set(ax, 'linewidth', 1.5);
                %                 xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
                %                 ylabel('\fontname {times new roman} Wavelength (pixel)', 'fontsize', 16);
            end
            cf = figure();
            ax = subplot(1, 1, 1);
            plot(theta, angular_1D);
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Sum of the magnitude response (arb.)', 'fontsize', 16);
        end
        
        function show_orientation_by_ID_RT_oneradius(obj, r_idx, radius, theta)
            % demonstrate the orientation extracted by Information diagram
            % of RT
            % use the inam image calcualated before
            % r_idx: the index of the radius selected
            % radius, theta: the same ratio for the multi-resolution RT
            % ***
            [lx, ly] = size(obj.C_scan_inam);
            % remove the edges
            image_orientation = NaN(lx, ly);
            for i = 1: lx
                for j = 1: ly
                    ID = squeeze(obj.RT_ID(i, j, r_idx, :));
                    %                     angular_1D = sum(ID, 1);
                    %                     [~, indAng] = max(angular_1D);
                    [~, indAng] = max(ID);
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
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \theta (\circ)');
        end
        
        function show_orientation_by_ID_RT(obj, radius, theta, imagename, angle_compens)
            % demonstrate the orientation extracted by Information diagram
            % of RT
            % use the inam image calcualated before
            % radius, theta: the same ratio for the multi-resolution RT
            % angle_compens: the angle used for composation
            % ***
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
            if ~isempty(obj.rear_mask)
                image_orientation(~obj.rear_mask) = NaN;
            end
            %
            if size(ID, 2) > 1
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
%             image_orientation = imgaussfilt(image_orientation, 1);
            %
            cf = figure('Name', ['orientation_image_RT_' num2str(r_max)]);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            shading flat;
            colormap hsv;
            set(gca, 'YDir', 'reverse');
            caxis([-90 90]);
            h = colorbar;
            h.Location = 'northoutside';
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \theta (\circ)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} x (pixel)', 'fontsize', 16);
            ylabel('\fontname {times new roman} y (pixel)', 'fontsize', 16);
            %
            figure('Name', ['radius_image_RT_']);
            ax = subplot(1, 1, 1);
            X = (1: size(image_radius, 2)); % pixels
            Y = (1: size(image_radius, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_radius);
            shading flat;
            colormap jet;
            set(gca, 'YDir', 'reverse');
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Radius (pixel)');
        end
        
        function show_inplane_direction_3D(obj, xslice, yslice, zslice, angle_compens)
            % demonstrate the 3d inplane_direction extracted by Gabor filter
            % demonstrate the ID_sum_overall
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % angle_compens: the compensation for angle
            % ***
            inph_ex = obj.Inplane_direction_3D;
            inph_ex = mod(inph_ex + angle_compens, 180) - 90;
            if isempty(inph_ex)
                inph_ex = obj.Inplane_direction_3D_ID;
            end
            y         = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x         = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            % ******************** show ply tracks
            cf = figure('Name', ['3d_orientation_RT_', 'R_' num2str(obj.radius_RT(1)), 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            x_idx = xslice * obj.fx / 1e3;
            y_idx = yslice * obj.fy / 1e3;
            mask  = (obj.col==x_idx) | (obj.row==y_idx);
            x_row = obj.row(mask) / obj.fx * 1e3;
            x_col = obj.col(mask) / obj.fy * 1e3;  % / obj.fy * 1e3 ;
            x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
            %             scatter3(ax, x_col, x_row, x_dep, ...
            %                 3, [122 122 121]/255, 'filled', ...
            %                 'DisplayName','Interply track');
            %             hold on;
            if sum(isnan(obj.front_I))~=0
                scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                    y, obj.front_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
                    3, 'black', 'filled', ...
                    'DisplayName','Front surface');
                hold on;
                scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                    y, obj.rear_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
                    3,  [0.75 0.75 0.75], 'filled', ...
                    'DisplayName','Rear surface');
                hold on;
                % select the yslice
                scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                    obj.front_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
                    3, 'black', 'HandleVisibility','off');
                hold on;
                scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                    obj.rear_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
                    3, [0.75 0.75 0.75], 'HandleVisibility','off');
                hold on;
            end
            % ******************* in-plane angle display *************
            % the 1st dimension is length(radius)
            if size(size(inph_ex), 1)==4 || size(size(inph_ex), 2)==4  % check dims
                z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
                [X, Y, Z] = meshgrid(x, y, z);
                % inph_visual(mask_phase_interp_visual==1) = NaN;
                for i = 1: size(inph_ex, 1)
                    h  = slice(ax, X, Y, Z, squeeze(inph_ex(i,:, :, :)), xslice, yslice, zslice);
                    hold on;
                    set(h, 'EdgeColor', 'none');
                    colormap hsv;
                    h  = colorbar;
                    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
                    set(ax, 'fontsize', 16);
                    set(ax, 'linewidth', 1.5);
                    xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
                    ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
                    zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
                    set(gca, 'ZDir', 'reverse');
                end
            else
                z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
                [X, Y, Z] = meshgrid(x, y, z);
                h         = slice(ax, X, Y, Z, inph_ex , xslice, yslice, zslice);
                hold on;
                set(h, 'EdgeColor', 'none');
                colormap hsv;
                h         = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \theta (\circ)');
                set(ax, 'fontsize', 16);
                set(ax, 'linewidth', 1.5);
                xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
                ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
                zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
                set(gca, 'ZDir', 'reverse');
            end
            %             lgd = legend('Interply', 'Front surface', 'Rear surface');
            lgd = legend('Front surface', 'Rear surface');
            xlim([x(1) x(end)]);
            ylim([y(1) y(end)]);
%             zlim([1 5.5])
            view([15 65 40]);
            % ************ Inplane_orientation_3D_overall ******
            inph_ex = obj.Inplane_direction_3D;
            if isempty(inph_ex)
                inph_ex     = obj.Inplane_direction_3D_overall_ID;
            end
            y           = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x           = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z           = (0: size(inph_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
            [X, Y, Z]   = meshgrid(x, y, z);
            inph_visual = inph_ex;
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf          = figure('Name', ['3d_orientation_RT_overall_' num2str(obj.radius_RT) '_' num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax          = subplot(1, 1, 1);
            h           = slice(ax, X, Y, Z, inph_visual, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap hsv;
            h           = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            % ************ Inplane_orientation_3D_overall_layerbylayer ******
            inph_ex     = obj.Inplane_direction_3D_overall_layer;
            y           = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x           = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z           = (0: size(inph_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
            [X, Y, Z]   = meshgrid(x, y, z);
            inph_visual = inph_ex;
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf          = figure('Name', ['3d_orientation_RT_overall_' num2str(obj.radius_RT) '_' num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax          = subplot(1, 1, 1);
            h           = slice(ax, X, Y, Z, inph_visual, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap hsv;
            h           = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
        end
        
        function show_inplane_direction_3D_meanstd(obj, xslice, yslice, zslice, plies, ratios, PropertyName, wavelength)
            % demonstrate the 3d inplane_direction extracted by Gabor filter
            % demonstrate the ID_sum_overall
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % ************ mean ************
            inph_ex                         = obj.Inplane_direction_3D;
            [lx, ly, ~]                     = size(inph_ex);
            Inplane_d_plywise_mean          = NaN(size(inph_ex));
            Inplane_d_plywise_std           = NaN(size(inph_ex));
            Inplane_d_plywise_mean_wholeply = NaN(size(inph_ex));
            wl_max                          = max(wavelength);
            plywise_direction               = NaN(size(inph_ex, 1), size(inph_ex, 2), plies, length(ratios));
            for p = 1:plies
                tic;
                % check the mean and std values of the orientations in one ply
                [~, ~, upper_index_bound]    = obj.define_plywise_inamCscan ...
                    (p, 0, PropertyName);
                [~, ~, lower_index_bound]    = obj.define_plywise_inamCscan ...
                    (p, 1, PropertyName);
                for ratio_index = 1:length(ratios)
                    for ii = round(wl_max / 2) + 1: round(lx - wl_max / 2)
                        for jj = round(wl_max / 2) + 1: round(ly - wl_max / 2)
                            plywise_direction(ii, jj, p, ratio_index) = mode(inph_ex(ii, jj, upper_index_bound(ii, jj): ...
                                lower_index_bound(ii, jj)), 'all');
                            Inplane_d_plywise_mean(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj)) ...
                                = mode(inph_ex(ii, jj, upper_index_bound(ii, jj): ...
                                lower_index_bound(ii, jj)), 'all');
                            Inplane_d_plywise_std(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj)) ...
                                = sum( abs(inph_ex(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj)) - ...
                                Inplane_d_plywise_mean(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj))), 'all', 'omitnan');
                        end
                    end
                end
                % assign the orientation for a ply
                for ii = round(wl_max / 2) + 1: round(lx - wl_max / 2)
                    for jj = round(wl_max / 2) + 1: round(ly - wl_max / 2)
                        Inplane_d_plywise_mean_wholeply(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj)) ...
                            = mode(plywise_direction(:, :, p, :), 'all');
                    end
                end
                toc;
            end
            % ********* display mean **************
            y         = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x         = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            cf        = figure('Name', ['3d_orientation_RT_mean_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax        = subplot(1, 1, 1);
            h         = slice(ax, X, Y, Z, Inplane_d_plywise_mean, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap hsv;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mean angle (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            % ********* std **************
            y         = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x         = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            cf        = figure('Name', ['3d_orientation_RT_mean_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax        = subplot(1, 1, 1);
            h         = slice(ax, X, Y, Z, Inplane_d_plywise_std , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} STD (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            % ********* mean whole ply **************
            y         = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x         = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            cf        = figure('Name', ['3d_orientation_RT_mean_wholeply_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax        = subplot(1, 1, 1);
            h         = slice(ax, X, Y, Z, Inplane_d_plywise_mean_wholeply , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} STD (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
        end
        
        function obj = extract_local_orientation_RT_3D_parallel(obj, PropertyName, radius, theta, ratios, sigma_denoise)
            % extract the local orientation image in the plane parallel to
            % surfaces
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % radius, theta: the ratio for the multi-resolution RT
            % ratios: the ratio determing the index by the front and back surfaces , unit: arb.
            % sigma_denoise: threshould for curvelet denoising, 0 means no application of denoising
            % ***
            % assign the info for saving
            obj.radius_RT             = radius;
            obj.theta_RT              = theta;
            obj.PropertyName_RT       = PropertyName;
            % start
            [lxc, lyc]                = size(obj.front_I);
            Inplane_direction         = NaN(lxc, lyc, size(obj.(PropertyName), 3));
            Inplane_direction_overall = NaN(size(obj.(PropertyName)));
            index_container           = NaN(lxc, lyc);
            r_max                     = max(radius);
            %             ID_sum_all = NaN(length(distances_to_front), length(radius), length(theta));
            for i = 1:length(ratios)
                ratio            = ratios(i); % the ratio to determine the index
                if i~=length(ratios)
                    ratio_next       = ratios(i+1);
                end
                index_container_prev = index_container;
                [~, C_scan_inam_para, index_container] = obj.define_parallel_inamCscan(ratio, PropertyName);
                [~, ~, index_container_next]           = obj.define_parallel_inamCscan(ratio_next, PropertyName);
                if sigma_denoise==0
                    C_scan_inam_para_denoise = C_scan_inam_para;
                else
                    C_scan_inam_para_denoise = fx_curvelet_denoise_enhanced(sigma_denoise, C_scan_inam_para);
                end
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
                            if maxval <= mean(anguler_1D) + std(anguler_1D) % there is or is not significant fibrous content.
                                continue;
                            end
                            if i==1
                                upper_index_bound = index_container(ii, jj);
                                lower_index_bound = round(index_container(ii, jj)/2 ...
                                    +index_container_next(ii, jj)/2);
                            elseif i==length(ratios)
                                upper_index_bound = round(index_container(ii, jj)/2 ...
                                    + index_container_prev(ii, jj)/2);
                                lower_index_bound = index_container(ii, jj);
                            else
                                upper_index_bound = round(index_container(ii, jj)/2 ...
                                    + index_container_prev(ii, jj)/2);
                                lower_index_bound = round(index_container(ii, jj)/2 ...
                                    +index_container_next(ii, jj)/2);
                            end
                            Inplane_direction(ii, jj, upper_index_bound: ...
                                lower_index_bound) = theta(I);
                        end
                    end
                end
                %                 disp('finished the RT:');
                %                 toc;
                % search for the max
                %                 ID_sum = zeros(size(ID, 3), size(ID, 4));
                %                 for ii = 1: lxc
                %                     for jj = 1: lyc
                %                         ID_pixel = squeeze(ID(ii, jj, :, :));
                %                         %                         front_index = obj.front_I(ii, jj);
                %                         %                         % special transverse when one dimension is 1;
                %                         %                         if size(ID_pixel, 2)==1
                %                         %                             ID_pixel = ID_pixel';
                %                         %                         end
                %                         %                         [~, I] = max(ID_pixel(:));
                %                         %                         [~, indAng] = ind2sub(size(ID_pixel),I);
                %                         %                         if front_index + dis_to_front < obj.rear_I(ii, jj) && front_index > 0
                %                         %                             Inplane_direction(ii, jj, round(front_index + dis_to_front - spacing_dis / 2)...
                %                         %                                 :round(front_index + dis_to_front + spacing_dis / 2)) = theta(indAng);
                %                         %                         end
                %                         % sum up the ID
                %                         ID_sum = ID_sum + ID_pixel;
                %                     end
                %                 end
                % end loop in one Cscan
                disp([num2str(ratios(i)), '/', num2str(ratios(end))]);
            end
            obj.Inplane_direction_3D         = Inplane_direction;
            obj.Inplane_direction_3D_overall = Inplane_direction_overall;
%             obj.ID_sum_overall               = ID_sum_all;
        end
          
        function obj = extract_local_orientation_RT_3D_plywise(obj, PropertyName, radius, theta, ratios, plies, sigma)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % radius, theta: the ratio for the multi-resolution RT
            % ratios: the ratios to define the plane in one ply
            % plies: number of the plies
            % sigma: the threshold for denoising 
            % ***
            % assign the info for saving
            obj.radius_RT             = radius;
            obj.theta_RT              = theta;
            obj.PropertyName_RT       = PropertyName;
            % start
            Inplane_direction         = NaN(length(radius), ...
                size(obj.(PropertyName), 1),...
                size(obj.(PropertyName), 2),...
                size(obj.(PropertyName), 3));
            Inplane_direction_overall = NaN(size(obj.(PropertyName)));
            index_container           = NaN(size(obj.front_I));
            ID_sum_all                = NaN(length(ratios)*plies, length(radius), length(theta));
            for p = 1:plies
                tic;
                for ratio_index = 1:length(ratios)
                    ratio                                 = ratios(ratio_index);
                    index_container_prev                  = index_container;
                    [~, C_scan_inam_ply, index_container] = obj.define_plywise_inamCscan(p, ratio, PropertyName);
                    if ratio_index==1
                        [~, ~, upper_index_bound]    = obj.define_plywise_inamCscan ...
                            (p, 0, PropertyName);
                        [~, ~, index_container_next] = obj.define_plywise_inamCscan ...
                            (p, ratios(ratio_index+1), PropertyName);
                        lower_index_bound            = round(index_container/2 ...
                            +index_container_next/2);
                    elseif ratio_index==length(ratios)
                        upper_index_bound            = round(index_container/2 ...
                            + index_container_prev/2);
                        [~, ~, lower_index_bound]    = obj.define_plywise_inamCscan ...
                            (p, 1, PropertyName);
                    else
                        [~, ~, index_container_next] = obj.define_plywise_inamCscan ...
                            (p, ratios(ratio_index+1), PropertyName);
                        upper_index_bound            = round(index_container/2 ...
                            + index_container_prev/2);
                        lower_index_bound            = round(index_container/2 ...
                            +index_container_next/2);
                    end
                    if sigma~=0
                        C_scan_inam_ply = fx_curvelet_denoise_enhanced(sigma, C_scan_inam_ply);
                    end
                    % Radon transform
%                     fill nan
                    C_scan_inam_ply    = fx_inpaint_nans(C_scan_inam_ply, 5);
                    % pad the C_scan_inam matrix
                    C_scan_inam_RT_pad = C_scan_inam_ply;
                    r_max              = max(radius);
                    [lxc, lyc]         = size(C_scan_inam_RT_pad);
                    ID                 = NaN(lxc, lyc, length(radius), length(theta));
                    % check the nan and adjust the radii
                    NaN_part           = isnan(C_scan_inam_RT_pad);
                    [r_nan, l_nan]     = find(NaN_part);
                    for ii = r_max + 1: lxc - r_max
                        for jj = r_max + 1: lyc - r_max
                            if obj.rear_I(ii, jj) <= index_container(ii, jj)% check rear surface
                                continue;
                            end
                            for r_idx = 1: length(radius)
                                r      = radius(r_idx);
                                % check the nan and adjust the radii
                                distances = sqrt((r_nan-ii).^2 + (l_nan-jj).^2);
                                if r > min(distances)
%                                     break;
                                    r = floor(min(distances));
                                end
                                if r < 4
                                    break;
                                end
                                center                = [jj, ii];
                                % creat a round mask parameters
                                [~, anguler_1D, ~, ~] = fx_Radonto1Dangular_correct(C_scan_inam_RT_pad, center, r, theta);
                                ID(ii, jj, r_idx, :)  = anguler_1D;
                                % search for the max
                                [~, I]                = max(anguler_1D);
                                if ~isnan(upper_index_bound(ii, jj)) && ...
                                        ~isnan(lower_index_bound(ii, jj))
                                    Inplane_direction(r_idx, ii, jj, upper_index_bound(ii, jj): ...
                                        lower_index_bound(ii, jj)) = theta(I);
                                end
                            end
                        end
                    end
                    % search for the max
%                     ID_sum = zeros(size(ID, 3), size(ID, 4));
%                     for ii = 1: lxc
%                         for jj = 1: lyc
%                             obj.ID_pixel = squeeze(ID(ii, jj, :, :));
%                             ID_sum = ID_sum + obj.ID_pixel;
%                         end
%                     end
                end
                toc;
                disp(['ply:', num2str(p), ' / ', num2str(plies(end))]);
            end
            obj.Inplane_direction_3D         = squeeze(Inplane_direction);
            obj.Inplane_direction_3D_overall = squeeze(Inplane_direction_overall);
            obj.ID_sum_overall               = ID_sum_all;
        end
 
        % ************************** curvelet denoise & curvelt orientaion extraction *****************
        function obj = curvelet_denoise_inamCscan(obj, sigma, imagename)
            % sigma: noise level, also a parameter to define the threshold
            % imagename: name of the image property
            % basis curvelet denoising from 'fdct_wrapping_demo_denoise'
            % fdct_wrapping_demo_basic.m -- Displays a curvelet both in the spatial and frequency domains.
            noisy_img = obj.(imagename);
            NaN_part  = isnan(noisy_img);
            X         = (1: size(noisy_img, 2)) / obj.fx * 1e3; % pixels
            Y         = (1: size(noisy_img, 1)) / obj.fy * 1e3;
            %             se                          = strel('disk',5);
            %             closeNaN_part               = imopen(NaN_part, se); % Morphologically close imag
            restored_img           = fx_curvelet_denoise_enhanced(sigma, noisy_img);
            restored_img(NaN_part) = NaN;
%             title('Original Image')
            cf      = figure('Name', 'x-y slice-plane of instantaneous amplitude'); 
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            clf; ax = subplot(1,1,1);
            pcolor(ax, X, Y, noisy_img);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap gray;
            h          = colorbar;
            h.Location = 'northoutside';
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            %
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
%             figure; clf; imagesc(restored_img); colormap gray; axis('image');
%             title('Restored image');
            obj.C_scan_inam_denoise = restored_img;
            colorbar;
        end
        
        function obj = compute_curvelet(obj, imagename)
            % imagename: name of the image property
            % Forward curvelet transform
            noisy_img                   = obj.(imagename);
            disp('Take curvelet transform: fdct_usfft');
            noisy_img(isnan(noisy_img)) = mean(noisy_img, [1 2], 'omitnan');
            tic; C     = fdct_wrapping(double(noisy_img), 0, 1); toc;
            % upsample the Coefficients by duplication
            % save as Information diagram
            [img_m, img_n] = size(noisy_img);
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
            obj.Curvelet_ID   = C_temp;
            obj.Curvelet_Coef = C;
        end
        
        function show_orientation_by_curveletID(obj, imagename)
            % demonstrate the orientation extracted by curvelet ID
            % use the 4D images calcualated before
            % imagename: name of the image property
            % orientation: the orientation for Gabor filter, to
            % ***
            % kernel_size_half =
            [lx, ly]          = size(obj.(imagename));
            image_orientation = NaN(lx, ly);
            image_scale       = NaN(lx, ly);
            Curvelet_ID_4d    = obj.Curvelet_ID;
            for i = 1: lx
                for j = 1: ly
                    ID                      = squeeze(Curvelet_ID_4d(i, j, 1:end, :));
                    [~, I]                  = max(ID(:));
                    [indwl, indAng]         = ind2sub(size(ID), I);
                    image_orientation(i, j) = -mod((45 + indAng/size(Curvelet_ID_4d, 4)*360), -180); % start from 45 degree
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
        
        % not yet correct
        function [C_scan_inam_profileuwp, obj]= define_Cscan_depthprofile_uwphase(obj, dis_to_front, PropertyName)
            % define the Cscan of inam by unwrapped phase in 3d
            % save the C_scan_inam
            % dis_to_front: the distance to the front surface , unit: us
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % ***
            if strcmp(PropertyName, 'img')
                inam = obj.(PropertyName);
            else
                inam = abs(obj.(PropertyName));
            end
            dis_to_front = dis_to_front * obj.fs / 1e6; % us to sample points
            C_scan_uwphase = NaN([size(inam, 1), size(inam, 2)]);
            inph_unwp = obj.inph_unwphase;
            fronts_index = obj.front_I;
            rears_index = obj.rear_I;
            % creat a mask for rear surface
            obj.rear_mask = obj.front_I + dis_to_front < obj.rear_I;
            for i = 1: size(inam, 1)
                for j = 1:size(inam, 2)
                    front_index = fronts_index(i, j);
                    if front_index + dis_to_front < rears_index(i, j) - 50
                        % here, set a threshold to define the rear surface.
                        C_scan_uwphase(i, j) = inph_unwp(i, j, round(front_index + dis_to_front));
                    end
                end
            end
            % calculate the average unwrapped phase
            C_scan_uwphase_ave = mean(C_scan_uwphase, 'all', 'omitnan');
            C_scan_inam_profileuwp = NaN([size(inam, 1), size(inam, 2)]);
            for idx = 1: size(inam, 1)
                for idy = 1:size(inam, 2)
                    if rears_index(idx, idy) <=0 || rears_index(idx, idy) > size(inph_unwp, 3) ...
                            || fronts_index(idx, idy) <=0 || fronts_index(idx, idy) > rears_index(idx, idy)
                        continue;
                    end
                    % track the phase by unwraped inst. phase
                    A_inph_unw = squeeze(inph_unwp(idx, idy, fronts_index(idx, idy):rears_index(idx, idy)));
                    % remove duplicates using 'unique' function.
                    [A_inph_unw, ~] = unique(A_inph_unw);
                    depindex = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                        C_scan_uwphase_ave) + fronts_index(idx, idy);
                    if depindex > 0 && depindex < rears_index(idx, idy)
                        C_scan_inam_profileuwp(idx, idy) = inam(idx, idy, round(depindex));
                    end
                end
            end
        end

        function obj = define_followUnwrapP_inamCscan(obj, dis_to_front, PropertyName)
            % demonstrate the profile of one inter-ply in 3d
            % save the C_scan_inam
            % dis_to_front: the distance to the front surface , unit: us
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % ***
            if strcmp(PropertyName, 'img')
                inam = obj.(PropertyName);
            else
                inam = abs(obj.(PropertyName));
            end
            dis_to_front = dis_to_front * obj.fs / 1e6; % us to sample points
            obj.C_scan_inam = NaN([size(inam, 1), size(inam, 2)]);
            C_scan_uwphase = NaN([size(inam, 1), size(inam, 2)]);
            C_scan_inam_parallel = NaN([size(inam, 1), size(inam, 2)]);
            inph_unwp = obj.inph_unwphase;
            fronts_index = obj.front_I;
            rears_index = obj.rear_I;
            % creat a mask for rear surface
            obj.rear_mask = obj.front_I + dis_to_front < obj.rear_I;
            for i = 1: size(inam, 1)
                for j = 1:size(inam, 2)
                    front_index = fronts_index(i, j);
                    if front_index + dis_to_front < rears_index(i, j) - 50
                        % here, set a threshold to define the rear surface.
                        C_scan_uwphase(i, j) = inph_unwp(i, j, round(front_index + dis_to_front));
                        C_scan_inam_parallel(i, j) = inam(i, j, round(front_index + dis_to_front));
                    end
                end
            end
            disp(dis_to_front);
            % calculate the average unwrapped phase
            C_scan_uwphase_ave = mean(C_scan_uwphase, 'all', 'omitnan');
            C_scan_profile_uwpave = NaN([size(inam, 1), size(inam, 2)]);
            C_scan_inam_profileuwp = NaN([size(inam, 1), size(inam, 2)]);
            for idx = 1: size(inam, 1)
                for idy = 1:size(inam, 2)
                    if rears_index(idx, idy) <=0 || rears_index(idx, idy) > size(inph_unwp, 3) ...
                            || fronts_index(idx, idy) <=0 || fronts_index(idx, idy) > rears_index(idx, idy)
                        continue;
                    end
                    % track the phase by unwraped inst. phase
                    A_inph_unw = squeeze(inph_unwp(idx, idy, fronts_index(idx, idy):rears_index(idx, idy)));
                    % remove duplicates using 'unique' function.
                    [A_inph_unw, ~] = unique(A_inph_unw);
                    depindex = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                        C_scan_uwphase_ave) + fronts_index(idx, idy);
                    C_scan_profile_uwpave(idx, idy) = depindex + fronts_index(idx, idy);
                    if depindex > 0 && depindex < rears_index(idx, idy)
                        C_scan_inam_profileuwp(idx, idy) = inam(idx, idy, round(depindex));
                    end
                end
            end
            % save the image as obj's property
            obj.C_scan_inam = C_scan_inam_profileuwp;
            %
            X = (1: size(C_scan_uwphase, 2)) / obj.fx * 1e3;
            Y = (1: size(C_scan_uwphase, 1)) / obj.fy * 1e3;
            figure('Name', ['Cscan_inam_parallel_' , PropertyName, '_', num2str(dis_to_front)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, X, Y, C_scan_inam_parallel);
            hold on;
            h = colorbar;
            %
            figure('Name', ['Cscan_unwphase_parallel_' , PropertyName, '_', num2str(dis_to_front)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, X, Y, C_scan_uwphase);
            hold on;
            h = colorbar;
            %
            figure('Name', ['Cscan_profile_unwphase' , PropertyName, '_', num2str(dis_to_front)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, X, Y, C_scan_profile_uwpave);
            hold on;
            h = colorbar;
            %
            figure('Name', ['Cscan_inam_unwphase' , PropertyName, '_', num2str(dis_to_front)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, X, Y, C_scan_inam_profileuwp);
            hold on;
            h = colorbar;
        end
        
        function obj = circle_rotate(obj, l, w, center, degree)
            % rotate a circle in the image obj.C_scan_inam
            % radii: radius of the circle mask
            % center: [j, i], the centers of the circle mask
            % degree: the degree of rotation
            % l: Length
            % w: Width
            % *****************
            % creat a round mask
            radii = max(w, l);
            mask = fx_createCirclesMask(obj.C_scan_inam, center, radii);
            x1 = center(2) - radii;
            x2 = center(2) + radii;
            y1 = center(1) - radii;
            y2 = center(1) + radii;
            % use the mask to select a round region
            C_scan_inam_mask = obj.C_scan_inam .* mask;
            C_scan_inam_mask = C_scan_inam_mask(x1: x2,y1: y2);
            C_scan_inam_mask = imrotate(C_scan_inam_mask, degree, 'bicubic', 'crop');
            % creat a ellipse mask
            [X, Y] = meshgrid(1:size(C_scan_inam_mask, 1), 1:size(C_scan_inam_mask, 2));
            %make a meshgrid: use the size of your image instead
            X0 = size(C_scan_inam_mask, 1) / 2;
            Y0 = size(C_scan_inam_mask, 2) / 2;
            ellipse = ((X - X0) / l).^2+((Y -  Y0) / w).^2<=1;
            C_scan_inam_mask_rotate = C_scan_inam_mask .* ellipse;
            %
            for i = x1:x2
                for j = y1:y2
                    if C_scan_inam_mask_rotate(i - x1 + 1, j - y1 + 1)~=0
                        obj.C_scan_inam(i, j) = C_scan_inam_mask_rotate(i - x1 + 1, j - y1 + 1);
                    end
                end
            end
            figure('Name', ['Cscan_RT_parallel_rotated', '_', num2str(degree)]);
            ax = subplot(1, 1, 1);
            X = (1: size(obj.C_scan_inam, 2)) / obj.fx * 1e3;
            Y = (1: size(obj.C_scan_inam, 1)) / obj.fy * 1e3;
            imagesc(ax, X, Y, obj.C_scan_inam);
            hold on;
            h = colorbar;
            colormap jet;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} inst. amp. (arb,)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function obj = compute_orientation_by_RT(obj, radiis, theta)
            % extract the inplane orientation by Radon transform
            % radiis: the radius array for RT
            % theta: angle space
            % ***
            % fillna
            C_scan_inam_RT = obj.C_scan_inam;
            C_scan_inam_RT(isnan(C_scan_inam_RT)) = mean(C_scan_inam_RT, 'all', 'omitnan');
            %             C_scan_inam_RT(isnan(C_scan_inam_RT)) = mean(C_scan_inam_RT, 'all', 'omitnan');
            [lxc, lyc] = size(C_scan_inam_RT);
            ID = NaN(lxc, lyc, length(radiis), length(theta));
            %             ID_gyy = NaN(lxc, lyc, length(radiis), length(theta));
            % pad the C_scan_inam matrix
            r_max = max(radiis);
            C_scan_inam_RT_pad = padarray(C_scan_inam_RT, [r_max r_max], 'symmetric');
            %             C_scan_inam_RT_pad = padarray(C_scan_inam_RT,[r_max r_max],mean(C_scan_inam_RT, 'all', 'omitnan'));
            % timer
            tic;
            for i = 1: lxc
                for j = 1: lyc
                    for r_idx = 1: length(radiis)
                        r = radiis(r_idx);
                        center = [j + r_max, i + r_max];
                        % creat a round mask parameters
                        [~, anguler_1D, ~, ~] = fx_Radonto1Dangular_correct(C_scan_inam_RT_pad, center, r, theta);
                        ID(i, j, r_idx, :)    = anguler_1D;
                        %                         ID_gyy(i, j, r_idx, :) = anguler_1D_gyy;
                    end
                end
                if mod(i, 10)==0
                    disp([num2str(i), '\', num2str(lxc)]);
                end
            end
            timeElapsed = toc;
            disp(['Compute_MultiResolution_RT ', num2str(timeElapsed)]);
            obj.RT_ID = ID;
        end

        function obj = extract_local_orientation_modifiedRT_3D_parallel(obj, PropertyName, radius, theta, distances_to_front)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % radius, theta: the ratio for the multi-resolution RT
            % distances_to_front: the distance to the front surface , unit: us
            % ***
            % assign the info for saving
            obj.radius_RT = radius;
            obj.theta_RT = theta;
            obj.distances_to_front_RT = distances_to_front;
            obj.PropertyName_RT = PropertyName;
            % start
            Inplane_direction = NaN(length(radius), size(obj.(PropertyName), 1),...
                size(obj.(PropertyName), 2),...
                size(obj.(PropertyName), 3));
            Inplane_direction_overall = NaN(size(obj.(PropertyName)));
            Inplane_direction_multireso = NaN(size(obj.(PropertyName)));
            inam = abs(obj.(PropertyName));
            spacing_dis = (distances_to_front(end) - distances_to_front(1)) ...
                / (length(distances_to_front) - 1) * obj.fs / 1e6;
            ID_sum_all = NaN(length(distances_to_front), length(radius), length(theta));
            for i = 1:length(distances_to_front)
                %                 % time clock
                %                 tic;
                dis_to_front = distances_to_front(i) * obj.fs / 1e6; % us to sample points
                C_scan_inam_para = NaN([size(inam, 1), size(inam, 2)]);
                for ii = 1: size(inam, 1)
                    for jj = 1:size(inam, 2)
                        front_index = obj.front_I(ii, jj);
                        if front_index + dis_to_front < obj.rear_I(ii, jj) && front_index > 0
                            C_scan_inam_para(ii, jj) = inam(ii, jj, round(front_index + dis_to_front));
                        end
                    end
                end
                % Radon transform
                % fillna
                C_scan_inam_para(isnan(C_scan_inam_para)) = mean(C_scan_inam_para, 'all', 'omitnan');
                % pad the C_scan_inam matrix
                r_max = max(radius);
                [lxc, lyc] = size(C_scan_inam_para);
                ID = NaN(lxc, lyc, length(radius), length(theta));
                C_scan_inam_RT_pad = padarray(C_scan_inam_para,[r_max r_max],mean(C_scan_inam_para, 'all', 'omitnan'));
                %                 disp('defined the profile:');
                %                 toc;
                for ii = 1: lxc
                    for jj = 1: lyc
                        front_index = obj.front_I(ii, jj);
                        for r_idx = 1: length(radius)
                            r = radius(r_idx);
                            center = [jj + r_max, ii + r_max];
                            % creat a round mask parameters
                            [~, anguler_1D, ~, ~] = fx_Radonto1Dangular_correct(C_scan_inam_RT_pad, center, r, theta);
                            ID(ii, jj, r_idx, :) = anguler_1D;
                            % search for the max
                            [~, I] = max(anguler_1D);
                            if front_index + dis_to_front < obj.rear_I(ii, jj) && front_index > 0
                                Inplane_direction(r_idx, ii, jj, round(front_index + dis_to_front - spacing_dis / 2)...
                                    :round(front_index + dis_to_front + spacing_dis / 2)) = theta(I);
                            end
                        end
                    end
                end
                %                 disp('finished the RT:');
                %                 toc;
                % search for the max
                ID_sum = zeros(size(ID, 3), size(ID, 4));
                for ii = 1: lxc
                    for jj = 1: lyc
                        ID_pixel = squeeze(ID(ii, jj, :, :));
                        front_index = obj.front_I(ii, jj);
                        % special transverse when one dimension is 1;
                        if size(ID_pixel, 2)==1
                            ID_pixel = ID_pixel';
                        end
                        [~, I] = max(ID_pixel(:));
                        [~, indAng] = ind2sub(size(ID_pixel),I);
                        if front_index + dis_to_front < obj.rear_I(ii, jj) && front_index > 0
                            Inplane_direction_multireso(ii, jj, round(front_index + dis_to_front - spacing_dis / 2)...
                                :round(front_index + dis_to_front + spacing_dis / 2)) = theta(indAng);
                        end
                        %                         sum up the ID
                        ID_sum = ID_sum + ID_pixel;
                    end
                end
                for ii = 1: lxc
                    for jj = 1: lyc
                        front_index = obj.front_I(ii, jj);
                        [~, I] = max(ID_sum(:));
                        [~, indAng] = ind2sub(size(ID_sum),I);
                        % save the summary ID
                        if front_index + dis_to_front < obj.rear_I(ii, jj) && front_index > 0
                            Inplane_direction_overall(ii, jj, round(front_index + dis_to_front - spacing_dis / 2)...
                                :round(front_index + dis_to_front + spacing_dis / 2)) = theta(indAng);
                        end
                    end
                end
                ID_sum_all(i, :, :) = ID_sum;
                %                 disp('summed up the ID:');
                %                 toc;
                disp([num2str(distances_to_front(i)), '/', num2str(distances_to_front(end)), 'us']);
            end
            obj.Inplane_direction_3D_modifyRT           = Inplane_direction;
            obj.Inplane_direction_3D_multireso_modifyRT = Inplane_direction_multireso;
            obj.Inplane_direction_3D_overall_modifyRT   = Inplane_direction_overall;
            obj.ID_sum_overall_modifyRT                 = ID_sum_all;
        end
        
        function extract_local_ID_AlongDepth_RT(obj, PropertyName, radius, theta, x, y, distances_to_front)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % radius, theta: the ratio for the multi-resolution RT
            % x, y: the index of the location (x, y);
            % ***
            inam = abs(obj.(PropertyName));
            [lx, ly, ~] = size(inam);
            ID = NaN(length(distances_to_front), length(radius), length(theta));
            for r_idx = 1: length(radius)
                for i = 1:length(distances_to_front)
                    dis_to_front = distances_to_front(i) * obj.fs / 1e6; % us to sample points
                    C_scan_inam_para = NaN([lx, ly]);
                    for ii = 1: size(inam, 1)
                        for jj = 1:size(inam, 2)
                            front_index = obj.front_I(ii, jj);
                            if front_index + dis_to_front < obj.rear_I(ii, jj) && front_index > 0
                                C_scan_inam_para(ii, jj) = inam(ii, jj, round(front_index + dis_to_front));
                            end
                        end
                    end
                    % Radon transform
                    % fillna
                    C_scan_inam_para(isnan(C_scan_inam_para)) = mean(C_scan_inam_para, 'all', 'omitnan');
                    % pad the C_scan_inam matrix
                    r                  = radius(r_idx);
                    C_scan_inam_RT_pad = padarray(C_scan_inam_para, [r r], mean(C_scan_inam_para, 'all', 'omitnan'));
                    center = [y + r, x + r];
                    % creat a round mask parameters
                    [~, anguler_1D, ~, ~] = fx_Radonto1Dangular(C_scan_inam_RT_pad, center, r, theta);
                    ID(i, r_idx, :) = anguler_1D;
                end
                % display the ID along the depth for each r_idx
                cf = figure('Name', ['ID' '_' 'r' num2str(r) '_' 'x' num2str(x) '_' 'y' num2str(y)]);
                set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
                ax = subplot(1, 1, 1);
                imagesc(ax, squeeze(ID(:, r_idx, :)));
                colormap jet;
                h = colorbar;
                set(get(h, 'Title'), 'string', ...
                    '\fontname {times new roman}\fontsize {16} Mag. of angular distribution from RT (arb.)');
                set(ax, 'fontsize', 16);
                set(ax, 'linewidth', 1.5);
                xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
                ylabel('\fontname {times new roman} Depth (pixel)', 'fontsize', 16);
                set(gca, 'ZDir', 'reverse');
            end
        end
        
        function show_inplane_direction_3D_modifyRT(obj, xslice, yslice, zslice)
            % demonstrate the 3d inplane_direction extracted by modified RT
            % demonstrate the ID_sum_overall
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % ***
            inph_ex = obj.Inplane_direction_3D_modifyRT;
            % the 1st dimension is length(radius)
            y = (0: size(inph_ex, 2) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 3) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 4) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            for i = 1: size(inph_ex, 1)
                cf = figure('Name', ['3d_orientation_RTmodify_' 'r' num2str(obj.radius_RT(i)) '_' 'xslice', num2str(xslice(1))]);
                set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
                ax = subplot(1, 1, 1);
                h = slice(ax, X, Y, Z, squeeze(inph_ex(i,:, :, :)), xslice, yslice, zslice);
                hold on;
                set(h, 'EdgeColor', 'none');
                colormap hsv;
                h = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
                set(ax, 'fontsize', 16);
                set(ax, 'linewidth', 1.5);
                xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
                ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
                zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
                set(gca, 'ZDir', 'reverse');
            end
            % ************ Inplane_orientation_3D_multiresolution ******
            inph_ex = obj.Inplane_direction_3D_multireso_modifyRT;
            y = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            cf = figure('Name', ['3d_orientation_RT_overall_' num2str(obj.radius_RT) '_' num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            h = slice(ax, X, Y, Z, inph_ex, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap hsv;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            % ************ Inplane_orientation_3D_overall ******
            inph_ex = obj.Inplane_direction_3D_overall_modifyRT;
            y = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            cf = figure('Name', ['3d_orientation_RT_overall_' num2str(obj.radius_RT) '_' num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            h = slice(ax, X, Y, Z, inph_ex, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap hsv;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            % demonstrate the ID_sum_overall
            % normalize each layer
            ID_sum_all = squeeze(obj.ID_sum_overall_modifyRT);
            for i = 1:size(ID_sum_all, 1)
                ID_sum_all(i, :) = ID_sum_all(i, :) / max(ID_sum_all(i, :));
            end
            cf = figure('Name', ['3d_ID_sum_overall' '_' num2str(obj.radius_RT)]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            imagesc(ax, squeeze(sum(ID_sum_all, 2)));
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', ...
                '\fontname {times new roman}\fontsize {16} Mag. of angular distribution from RT (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Depth (pixel)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
        end
        
        % EDA for RT
        function show_one_RT_2dfft(obj, layer, ratio, PropertyName, center, radii)
            % demonstrate the profile of one inter-ply in 3d
            % layer: number of the layer
            % ratio: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % ***
            % the ratio should be 0 ~ 1
            if ratio < 0 || ratio > 1
                msg = "The 'raio' should be between 0 and 1";
                error(msg);
            end
            %             centers = [j, i];
            %             z_index = front_ave;
            %             while z_index<rear_ave
            %                 inam_C_scan = squeeze(inam(:, :, z_index));
            %                 % creat a round mask parameters
            %                 [~, ~, max_angle_I, ~, theta] = fx_Radonto1Dangular(inam_C_scan, centers, radii);
            inam = abs(obj.(PropertyName));
            inph_ex = obj.mask_plytrack;
            profile_layer = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_inam = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    inph_ex(i, j, [round(obj.front_I(i, j)) round(obj.rear_I(i, j))]) = 1; % add the front_surface and the rear_surface
                    k = find(inph_ex(i, j, :), layer + 1);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                    % extract two inter-plies besides the ply
                    if length(k) > layer
                        k1 = k(layer);
                        k2 = k(layer + 1);
                        % use the ratio to calculate the position in the ply
                        z_index = round(k1 + (k2 - k1) * ratio);
                        C_scan_inam(i, j) = inam(i, j, z_index);
                        profile_layer(i, j) = z_index / obj.fs * 1e6;
                    else % set the last inter-ply as NaN otherwise
                        profile_layer(i, j) = NaN;
                    end
                end
            end
            % creat a round mask
            mask = fx_createCirclesMask(C_scan_inam, center, radii);
            % fillna
            % C_scan_inam =  fx_inpaint_nans(C_scan_inam, 5);
            x1 = center(2) - radii;
            x2 = center(2) + radii;
            y1 = center(1) - radii;
            y2 = center(1) + radii;
            % use the mask to select a round region
            C_scan_inam_mask = C_scan_inam .* mask;
            C_scan_inam_mask = C_scan_inam_mask(x1: x2,y1: y2);
            tic;
            % now make 2D fft of original image
            fft2D = fft2(C_scan_inam_mask);
            fft2D_shifted = fftshift(fft2D);
            % Gaussian Filter Response Calculation
            [M, N]=size(C_scan_inam_mask); % image size
            R = 2^0; % filter size parameter
            X = 0: M - 1;
            Y = 0: N - 1;
            [X, Y]=meshgrid(Y, X);
            Cx = 0.5 * N;
            Cy = 0.5 * M;
            Lo_h = exp(-((X - Cx).^2 + (Y - Cy).^2)./(2 * R).^2);
            Hi = 1 - Lo_h; % High pass filter=1-low pass filter
            R_l = 2^12;
            Lo = exp(-((X - Cx).^2 + (Y - Cy).^2)./(2 * R_l).^2);
            fft2D_shifted_filter = abs(fft2D_shifted .* Hi .* Lo);
            %
            timeElapsed = toc;
            disp(['2d fft: ', num2str(timeElapsed)]);
            figure('Name', ['FFT2D_amp_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, log(1+abs(fft2D_shifted_filter)) );
            figure('Name', ['FFT2D_phase_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, angle(fft2D_shifted_filter));
            %             % transfer to polar coordinate
            %             x = -radii: radii;
            %             y = -radii: radii;
            %             [theta, rho] = cart2pol(x, y);
            %             figure('Name', ['FFT2D_polar_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            %             ax = subplot(1, 1, 1);
            %             imagesc(ax, theta, rho, log(1+abs(fft2D_shifted)) );
        end
        
        function compare_one_RT_fitplytrack(obj, layer, ratio, PropertyName, center, radii)
            % demonstrate the profile of one inter-ply in 3d
            % layer: number of the layer
            % ratio: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % ***
            % the ratio should be 0 ~ 1
            if ratio < 0 || ratio > 1
                msg = "The 'raio' should be between 0 and 1";
                error(msg);
            end
            inam = abs(obj.(PropertyName));
            inph_ex = obj.mask_plytrack;
            profile_layer = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            C_scan_inam = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    inph_ex(i, j, [round(obj.front_I(i, j)) round(obj.rear_I(i, j))]) = 1; % add the front_surface and the rear_surface
                    k = find(inph_ex(i, j, :), layer + 1);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                    % extract two inter-plies besides the ply
                    if length(k) > layer
                        k1 = k(layer);
                        k2 = k(layer + 1);
                        % use the ratio to calculate the position in the ply
                        z_index = round(k1 + (k2 - k1) * ratio);
                        C_scan_inam(i, j) = inam(i, j, z_index);
                        profile_layer(i, j) = z_index / obj.fs * 1e6;
                    else % set the last inter-ply as NaN otherwise
                        profile_layer(i, j) = NaN;
                    end
                end
            end
            % Radon transform
            % creat a round mask parameters
            [R, anguler_1D, max_angle_I, xp, theta] = fx_Radonto1Dangular(C_scan_inam, center, radii);
            % plot the R-theta image
            figure('Name', ['Radon_transform_R_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, theta, xp, R);
            hold on;
            h = colorbar;
            % plot the 1D anguler (sum of gradient by central-different kernel)
            figure('Name', ['Radon_transform_1Danguler_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            plot(theta, anguler_1D, 'linewidth', 2);
            % find the parallel C-scan of the layer
            z_index = round(mean(profile_layer * obj.fs / 1e6, 'all', 'omitnan'));
            C_scan_inam_2 = inam(:, :, z_index);
            C_scan_inam_2(isnan(C_scan_inam)) = NaN;
            % creat a round mask parameters
            [R, anguler_1D, max_angle_I, xp, theta] = fx_Radonto1Dangular(C_scan_inam_2, center, radii);
            % plot the R-theta image
            figure('Name', ['Radon_transform_R_withoutplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, theta, xp, R);
            hold on;
            h = colorbar;
            % plot the 1D anguler (sum of gradient by central-different kernel)
            figure('Name', ['Radon_transform_1Danguler_withoutplytrack_' , PropertyName, '_', num2str(layer), '_', num2str(ratio), '_', num2str(center), '_', num2str(radii)]);
            plot(theta, anguler_1D, 'linewidth', 2);
        end
        
        function monogenic_signal_oneCscan_fitplytrack(obj, cw, center, radii)
            % demonstrate the monogenic signal of one inter-ply which has
            % been extracted in 'RT'
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % cw:  Centre-wavelengths in pixel units
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % ***
            % creat a round mask
            mask = fx_createCirclesMask(obj.C_scan_inam, center, radii);
            % fillna
            %  C_scan_inam =  fx_inpaint_nans(C_scan_inam, 5);
            x1 = center(2) - radii;
            x2 = center(2) + radii;
            y1 = center(1) - radii;
            y2 = center(1) + radii;
            % use the mask to select a round region
            C_scan_inam_mask = obj.C_scan_inam .* mask;
            C_scan_inam_mask = C_scan_inam_mask(x1: x2,y1: y2);
            I = fx_inpaint_nans(C_scan_inam_mask, 5);
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
            % Local phase (calculated on a per-scale basis)
            LP = localPhase(m1, m2, m3);
            % Local orientation (calculated on a per-scale basis)
            % Only need to pass the odd parts (m2,m3) as even part (m1) is irrelevant
            LO = localOrientation(m2, m3);
            % Display
            cf = figure();
            imagesc(I), axis image, axis off;
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            % colormap gray
            title('Test Image') ;
            figure();
            imagesc(reshape(LE, Y, [])), axis image, axis off, colormap gray;
            title('Local Energy Over Scales');
            figure();
            imagesc(reshape(LP, Y, [])), axis image, axis off, colormap gray;
            title('Local Phase Over Scales');
            figure();
            imagesc(reshape(LO, Y, [])), axis image, axis off; colorbar;
            title('Local Orientation Over Scales (radians)');
        end
        
        function RT_monogenic_signal_oneCscan_fitplytrack(obj, PropertyName, cw,  center, radii)
            % demonstrate the monogenic signal of one inter-ply
            % use the C_scan_inam defined before
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % cw:  Centre-wavelengths in pixel units
            C_scan_inam_mono = obj.C_scan_inam;
            mask_NaN = isnan(C_scan_inam_mono);
            % creat a round mask
            mask = fx_createCirclesMask(C_scan_inam_mono, center, radii);
            %  C_scan_inam =  fx_inpaint_nans(C_scan_inam, 5);
            x1 = center(2) - radii;
            x2 = center(2) + radii;
            y1 = center(1) - radii;
            y2 = center(1) + radii;
            C_scan_inam_mono(mask_NaN) = mean(C_scan_inam_mono, 'all', 'omitnan');
            [Y, X] = size(C_scan_inam_mono);
            %
            filtStruct = createMonogenicFilters(Y, X, cw, 'lg', 0.55);
            % Now we can use this structure to find the monogenic signal for the image
            [m1, m2, m3] = monogenicSignal(C_scan_inam_mono, filtStruct);
            %
            % Local energy (calculated on a per-scale basis)
            LE = localEnergy(m1, m2, m3);
            % Local phase (calculated on a per-scale basis)
            LP = localPhase(m1, m2, m3);
            % Local orientation (calculated on a per-scale basis)
            % Only need to pass the odd parts (m2,m3) as even part (m1) is irrelevant
            LO = localOrientation(m2, m3);
            % Radon transform
            % creat a round mask parameters
            LO(mask_NaN) = 0;
            LP(mask_NaN) = 0;
            [R, anguler_1D, max_angle_I, xp, theta] = fx_Radonto1Dangular(LO, center, radii);
            % plot the R-theta image
            figure('Name', ['Radon_transform_R_' , PropertyName, '_', num2str(center), '_', num2str(radii)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, theta, xp, R);
            hold on;
            h = colorbar;
            % plot the 1D anguler (sum of gradient by central-different kernel)
            figure('Name', ['Radon_transform_1Danguler_' , PropertyName, '_', num2str(center), '_', num2str(radii)]);
            plot(theta, anguler_1D, 'linewidth', 2);
            figure('Name', ['LO_RT_fitplytrack_' , PropertyName, '_']);
            ax = subplot(1, 1, 1);
            X = (1: size(LO, 1)) / obj.fx * 1e3;
            Y = (1: size(LO, 2)) / obj.fy * 1e3;
            imagesc(ax, X, Y, LO);
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function show_one_Radontransform_filter(obj, PropertyName, center, radii, R_lo_H, R_Lo)
            % demonstrate the RT of one inter-ply combined with filter
            % use the C_scan_inam defined before
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % R_lo_H: the filter size parameter for high pass fitler (1 - Lo_h)
            % R_Lo: the filter size parameter for high pass fitler (1 - Lo_h)
            % ***
            % fillna
            C_scan_inam_2DFFT = obj.C_scan_inam;
            C_scan_inam_2DFFT(isnan(C_scan_inam_2DFFT)) = mean(C_scan_inam_2DFFT, 'all', 'omitnan');
            % bandpass fitler in frequency domain
            % now make 2D fft of original image
            fft2D = fft2(C_scan_inam_2DFFT);
            fft2D_shifted = fftshift(fft2D);
            % Gaussian Filter Response Calculation
            [M, N] = size(C_scan_inam_2DFFT); % image size
            X = 0: M - 1;
            Y = 0: N - 1;
            [X, Y]=meshgrid(Y, X);
            Cx = 0.5 * N;
            Cy = 0.5 * M;
            Lo_h = exp(-((X - Cx).^2 + (Y - Cy).^2)./(2 * R_lo_H).^2);
            Hi = 1 - Lo_h; % High pass filter=1-low pass filter
            Lo = exp(-((X - Cx).^2 + (Y - Cy).^2)./(2 * R_Lo).^2);
            fft2D_shifted_filter = fft2D_shifted .* Hi .* Lo;
            fft2D_ishift = ifftshift(fft2D_shifted_filter);
            C_scan_inam_2DFFT = abs(ifft2(fft2D_ishift));
            % Radon transform
            % creat a round mask parameters
            [R, anguler_1D, max_angle_I, xp, theta] = fx_Radonto1Dangular(C_scan_inam_2DFFT, center, radii);
            % plot the R-theta image
            figure('Name', ['Radon_transform_R_' , PropertyName, '_', num2str(center), '_', num2str(radii)]);
            ax = subplot(1, 1, 1);
            imagesc(ax, theta, xp, R);
            hold on;
            h = colorbar;
            % plot the 1D anguler (sum of gradient by central-different kernel)
            figure('Name', ['Radon_transform_1Danguler_' , PropertyName, '_', num2str(center), '_', num2str(radii)]);
            plot(theta, anguler_1D, 'linewidth', 2);
            % creat a round mask
            mask = fx_createCirclesMask(C_scan_inam_2DFFT, center, radii);
            x1 = center(2) - radii;
            x2 = center(2) + radii;
            y1 = center(1) - radii;
            y2 = center(1) + radii;
            % use the mask to select a round region
            C_scan_inam_mask = C_scan_inam_2DFFT .* mask;
            C_scan_inam_mask = C_scan_inam_mask(x1: x2,y1: y2);
            figure('Name', ['Cscan_RT_fitplytrack_' , PropertyName]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_mask, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan_inam_mask, 2)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_mask);
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        % ************************* local In-plane direction by logGabor filter *********
        function obj = show_logGabor_filter(obj, PropertyName, wavelength, orientation)
            % demonstrate and save the logGabor filtered images
            % use the C_scan_inam defined before
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ***
            gaborArray = gabor(wavelength, orientation);
            % fillna
            C_scan_inam_LG = obj.C_scan_inam;
            C_scan_inam_LG(isnan(C_scan_inam_LG)) = mean(C_scan_inam_LG, 'all', 'omitnan');
            %
            figure('Name', ['Cscan_fitplytrack_' , PropertyName]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_LG, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan_inam_LG, 2)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_LG);
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % Apply filters to C-scan image.
            [obj.gaborMag, obj.gaborPha]= imgaborfilt(C_scan_inam_LG, gaborArray);
            % normalized the magnitude.
            for i = 1:length(gaborArray)
                BW = gaborArray(i).SpatialFrequencyBandwidth;
                sigmaX = gaborArray(i).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                sigmaY = sigmaX ./ gaborArray(i).SpatialAspectRatio;
                obj.gaborMag(:, :, i) = obj.gaborMag(:, :, i)  / (2 * sigmaX * sigmaY * pi);
            end
            obj.frow = length(wavelength);
            obj.fcol = length(orientation);
            cf = figure;
            for p = 1: obj.frow * obj.fcol
                subplot(obj.fcol, obj.frow, p);
                imshow(obj.gaborMag(:, :, p), []);
                theta = gaborArray(p).Orientation;
                lambda = gaborArray(p).Wavelength;
                title(sprintf('Orientation=%d, Wavelength=%d', theta, lambda));
                h = colorbar;
            end
            cf = figure;
            for p = 1: obj.frow * obj.fcol
                subplot(obj.fcol, obj.frow, p);
                imshow(obj.gaborPha(:, :, p), []);
                theta = gaborArray(p).Orientation;
                lambda = gaborArray(p).Wavelength;
                title(sprintf('Orientation=%d, Wavelength=%d', theta, lambda));
                h = colorbar;
            end
        end
        
        function obj = compute_logGabor_filter_withoutFig(obj, PropertyName, wavelength, orientation, SFB, SAR, imagename)
            % demonstrate and save the logGabor filtered images
            % use the C_scan_inam defined before
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % SFB, SAR: SpatialFrequencyBandwidth and SpatialAspectRatio
            % imagename: the name of the 2d image
            % ***
            gaborArray     = gabor(wavelength, orientation, ...
                'SpatialFrequencyBandwidth', SFB, 'SpatialAspectRatio', SAR);
            % find the amp. of the front surface
            %             Ascan_ave                         = squeeze(mean(obj.(PropertyName), [1 2], 'omitnan'));
            %             Ascan_ave_amp                     = abs(Ascan_ave);
            %             front_amp                         = Ascan_ave_amp(round(mean(obj.front_I, [1 2], 'omitnan')));
            %
            C_scan_inam_LG = obj.(imagename);
            NaN_part       = isnan(C_scan_inam_LG);
            % fill nan
            if sum(isnan(NaN_part))>0
                C_scan_inam_LG = fx_inpaint_nans(C_scan_inam_LG, 5);
            end
            % pad the C_scan_inam matrix
            %             wl_max                            = max(wavelength);
            C_scan_inam_RT_pad                = C_scan_inam_LG; % no padding
            %             C_scan_inam_RT_pad = padarray(C_scan_inam_LG, [round(wl_max / 2) round(wl_max / 2)], 'replicate');
            %             C_scan_inam_RT_pad                    = padarray(C_scan_inam_LG, [round(wl_max / 2) round(wl_max / 2)], mean(C_scan_inam_LG, 'all', 'omitnan'));
            %
            figure('Name', ['Cscan_fitplytrack_', PropertyName]);
            ax = subplot(1, 1, 1);
            X  = (1: size(C_scan_inam_RT_pad, 1)) / obj.fx * 1e3;
            Y  = (1: size(C_scan_inam_RT_pad, 2)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_RT_pad);
            hold on;
            h  = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % Apply filters to C-scan image. **replication padding**
            fprintf('Gabor filteration: start');
            [obj.gaborMag, obj.gaborPha] = imgaborfilt(C_scan_inam_RT_pad, gaborArray);
            % normalized the magnitude.
            for i = 1:length(gaborArray)
                BW                    = gaborArray(i).SpatialFrequencyBandwidth;
                sigmaX                = gaborArray(i).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                sigmaY                = sigmaX ./ gaborArray(i).SpatialAspectRatio;
                gaborMag_i            = squeeze(obj.gaborMag(:, :, i));
                gaborMag_i            = gaborMag_i  / (2 * sigmaX * sigmaY * pi);
                gaborMag_i(NaN_part)  = NaN;
                obj.gaborMag(:, :, i) = gaborMag_i;
                clc;
                fprintf('Gabor filteration: %0.2f%%\n',100*i/length(gaborArray));
            end
            obj.frow = length(wavelength);
            obj.fcol = length(orientation);
        end
        
        function show_ID_logGabor_filter(obj, x, y, wavelength, orientation)
            % demonstrate the information diagram of one point in the logGabor filtered images
            % use the logGabor filtered images calcualated before
            % x, y: index of the point in the image unit : pixel
            % wavelength, orientation: the same ratio for the logGabor filter
            % ***
            wl_max = max(wavelength);
            ID = squeeze(obj.gaborMag(x + round(wl_max / 2), y + round(wl_max / 2), :));
            ID = reshape(ID, [obj.frow, obj.fcol]);
            angular_1D = sum(ID, 1);
            % shift the angular_1D by 90 degrees (half length) in x axis
            %             angular_1D_shift = circshift(angular_1D, round(length(angular_1D)/2));
            %             angular_1D_sub = angular_1D - angular_1D_shift;
            [~, indAng] = max(angular_1D);
            disp(orientation(indAng));
            % fillna
            C_scan_inam_LG = obj.C_scan_inam;
            C_scan_inam_LG(isnan(C_scan_inam_LG)) = mean(C_scan_inam_LG, 'all', 'omitnan');
            %
            figure('Name', ['Cscan_fitplytrack_forlogGabor' ]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_LG, 2)) / obj.fx * 1e3; % default X = (1: size(C_scan_inam_LG, 2))  here
            Y = (1: size(C_scan_inam_LG, 1)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_LG);
            hold on;
            scatter(ax, y/ obj.fx * 1e3, x / obj.fy * 1e3, 30, 'r', 'filled'); % reversed x, y here
            colormap jet;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp.(arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % demonstrate the ID and 1D angular
            if length(wavelength) > 1
                cf = figure();
                ax = subplot(1, 1, 1);
                imagesc(orientation, wavelength, ID);
                h = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
                set(ax, 'fontsize', 16);
                set(ax, 'linewidth', 1.5);
                xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
                ylabel('\fontname {times new roman} Wavelength (pixel)', 'fontsize', 16);
            end
            %             cf = figure();
            %             ax = subplot(1, 1, 1);
            %             plot(orientation, angular_1D_sub);
            %             set(ax, 'fontsize', 16);
            %             set(ax, 'linewidth', 1.5);
            %             xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Sum of the magnitude response (arb.)', 'fontsize', 16);
        end
        
        function show_orientation_by_ID(obj, wavelength, orientation)
            % demonstrate the orientation extracted by Information diagram
            % use the logGabor filtered images calcualated before
            % orientation: the orientation for Gabor filter, to
            % determine the angle space
            % wavelength: the wavelength for Gabor filter, to determing the size of kernel for removing the edges
            % ***
            %             kernel_size_half =
            [lx, ly] = size(obj.C_scan_inam);
            % remove the edges
            wl_max = max(wavelength);
            image_orientation = NaN(lx, ly);
            for i = 1: lx
                for j = 1: ly
                    ID                      = squeeze(obj.gaborMag(i + round(wl_max / 2), j + round(wl_max / 2), :));
                    ID                      = reshape(ID, [obj.frow, obj.fcol]);
                    angular_1D              = sum(ID, 1);
                    % shift the angular_1D by 90 degrees (half length) in x axis
                    angular_1D_shift        = circshift(angular_1D, round(length(angular_1D)/2));
                    angular_1D_sub          = angular_1D - angular_1D_shift;
                    [~, indAng]             = max(angular_1D_sub);
                    image_orientation(i, j) = orientation(indAng);
                end
            end
            if ~isempty(obj.rear_mask)
                image_orientation = image_orientation .* obj.rear_mask;
            end
            %
            figure('Name', ['orientation_image_', num2str(wavelength(1))]);
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            imagesc(ax, X, Y, image_orientation);
            colormap(jet);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
        end
        
        function show_orientation_by_ID_allwl(obj, wavelength, orientation, K, imagename, angle_compens)
            % demonstrate the orientation extracted by Information diagram
            % use the logGabor filtered images calcualated before
            % orientation: the orientation for Gabor filter, to
            % determine the angle space
            % wavelength: the wavelength for Gabor filter, to determing the size of kernel for removing the edges
            % K: the smoothing applied to the Gabor magnitude responses.
            % imagename: the name of the 2d image property
            % angle_compens: angle for compensation
            % ***
            % kernel_size_half =
            [lx, ly]          = size(obj.(imagename));
            % remove the edges
            wl_max            = round(max(wavelength)/2);
%             wl_min            = round(min(wavelength)); % !! change this to adjust the edge cutting
            wl_min            = wl_max; % comment out to fit the image edges
            image_orientation = NaN(lx, ly);
            image_wavelength  = NaN(lx, ly);
            ID_sum            = zeros(length(wavelength), length(orientation));
            gabormag          = obj.gaborMag;
            g                 = gabor(wavelength,orientation);
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
            % check the nan and adjust the ID
            NaN_part       = isnan(gabormag);
            NaN_part       = sum(NaN_part, 3);
            [r_nan, l_nan] = find(NaN_part);
            % calculate the matrix of the distances to nan part
            [X, Y]         = meshgrid(1:ly, 1:lx);
            distances_pow2 = max(wavelength)^2 * ones(lx, ly);
            for i = 1:length(r_nan)
                distances_i_pow2 = (r_nan(i)-X).^2 + (l_nan(i)-Y).^2;
                distances_pow2   = min(distances_pow2, distances_i_pow2);
            end
            for i = round(wl_min)+1: lx - round(wl_min)-1
                for j = round(wl_min)+1: ly - round(wl_min)-1
                    ID                      = squeeze(gabormag(i, j, :));
                    if isnan(ID) % invalid part has been assigned NaN
                        image_orientation(i, j) = NaN;
                        image_wavelength(i, j)  = NaN;
                        continue;
                    end
                    ID        = reshape(ID, [obj.frow, obj.fcol]);
                    % check the nan and adjust the ID
                    %                     distances = sqrt((r_nan-i).^2 + (l_nan-j).^2);
                    %                     ID_bound  = fx_binarySearch(...
                    %                         wavelength, length(wavelength), min(distances));
                    %                     if ID_bound==1 || wavelength(ID_bound-1) <= 4
                    %                         continue;
                    %                     else
                    %                         ID = ID(1:ID_bound-1, :);
                    %                     end
                    % %wavelength is 2:1:max!!
                    % *********** comment out below for better performance on the ****edge !
                    %                     ID_bound = round(sqrt(distances_pow2(i, j)));
                    %                     if ID_bound <= 4
                    %                         continue;
                    %                     elseif ID_bound < wl_max*2 || ... % consider the edge of the whole image
                    %                             i < wl_max*2 || j < wl_max*2 || lx-i < wl_max*2 || ly-j < wl_max*2
                    %                         % consider the edge of the whole image
                    %                         ID_bound     = min([ID_bound, i, j, lx-i , ly-j ]);
                    %                         ID_bound_ind = fx_binarySearch(wavelength, length(wavelength), ID_bound)-1;
                    %                         ID           = ID(1:ID_bound_ind, :);
                    %                     end
                    %                     ID_shift = circshift(ID, round(size(ID, 2) / 2), 2);
                    % ****************** %
                    %                     ID_sub = ID - ID_shift;
                    [C, I]                  = max(ID(:));
                    if C < mean(ID, [1,2]) + std(ID, [], 'all')
                        continue;
                    end
                    [indwl, indAng]         = ind2sub(size(ID),I);
                    image_orientation(i, j) = mod(orientation(indAng)+angle_compens, 180) - 90;
                    image_wavelength(i, j)  = wavelength(indwl);
                    % sum up ID
                    if size(ID)==size(ID_sum)
                        ID_sum = ID_sum + ID;
                    end
                end
                clc;
                fprintf('ID imaging: %0.2f%%\n', 100*i /(lx - 2*round(wl_min)));
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
% %             display mean + std
%             disp([num2str(mean(image_orientation, 'all', 'omitnan')) ...
%                 ' + ' num2str(std(image_orientation, 0, 'all', 'omitnan'))]);
            %************************* Gaussian fitting ********************
            Idof_edge    = 0:180;
            Idof_N       = histcounts(image_orientation, Idof_edge);
            centers      = [0 45 90 135];
            sigmas       = [8 8 8 8];
            amplitudes   = [1 1 1 1];
            %             centers      = [45 90 135];
            %             sigmas       = [8 8 8];
            %             amplitudes   = [1 1 1];
            numGaussians = 4;
            [amp, parameter, ~] = fx_multi_Gaussianfitting(...
                numGaussians, Idof_edge(1:180), Idof_N, centers, sigmas, amplitudes);
            [~, amp_I]      = max(amp'.*parameter(2:2:end));
            m_fiber_angle   = parameter(amp_I*2-1);
            %             std_fiber_angle = amp(amp_I);
            %
            std_fiber_angle = std(image_orientation, 0, 'all','omitnan');
            disp(['mean:' num2str(m_fiber_angle) '. std:' num2str(std_fiber_angle)]);
            % ********** plot 1D distribution
            % normalize
            Idof_N = Idof_N / sum(Idof_N);
            % log
            figure('Name', ['angle_distribution_ID' '_' imagename]);
            set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
            h1 = bar(-89:90, Idof_N);
            xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
            %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
            set(gca, 'linewidth', 2);
            % %
            cf = figure('Name', ['orientation_image_', num2str(K), imagename, '_', num2str(angle_compens)]);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            X  = (1: size(image_orientation, 2)) / obj.fx * 1e3; % pixels
            Y  = (1: size(image_orientation, 1)) / obj.fy * 1e3; % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap hsv;
            caxis([-90 90]);
            h          = colorbar;
            h.Location = 'northoutside';
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \theta (\circ)');
            set(h,'YTick', [-90 -45 0 45 90]); % set ticks
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            xlim([min(X(:)) max(X(:))]);
            ylim([min(Y(:)) max(Y(:))]);
            %
            cf = figure('Name', ['wavelength_image', num2str(K), imagename, '_', num2str(angle_compens)]);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            X  = (1: size(image_wavelength, 2)) / obj.fx * 1e3; % pixels
            Y  = (1: size(image_wavelength, 1)) / obj.fy * 1e3; % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_wavelength);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap jet;
            h          = colorbar;
            h.Location = 'northoutside';
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} \lambda (pixels)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            xlim([min(X(:)) max(X(:))]);
            ylim([min(Y(:)) max(Y(:))]);
        end
        
        function image_orientation = show_orientation_by_ID_onewl(obj, wl_idx, wavelength, orientation)
            % demonstrate the orientation extracted by Information diagram
            % of RT
            % use the inam image calcualated before
            % r_idx: the index of the wavelength selected
            % orientation: the orientation for Gabor filter, to determine the angle space
            % wavelength: the wavelength for Gabor filter, to determing the size of kernel for removing the edges
            % ***
            [lx, ly] = size(obj.C_scan_inam);
            % remove the edges
            image_orientation = NaN(lx, ly);
            wl_max = max(wavelength);
            for i = 1: lx
                for j = 1: ly
                    ID = squeeze(obj.gaborMag(i + round(wl_max / 2), j + round(wl_max / 2), :));
                    ID = reshape(ID, [obj.frow, obj.fcol]);
                    % angular_1D = sum(ID, 1);
                    % [~, indAng] = max(angular_1D);
                    [~, indAng] = max(ID(wl_idx, :));
                    image_orientation(i, j) = orientation(indAng);
                end
            end
            disp(wavelength(wl_idx));
            if ~isempty(obj.rear_mask)
                image_orientation(~obj.rear_mask) = NaN;
            end
            %
            figure('Name', ['orientation_image_LG_' num2str(wavelength(wl_idx))]);
            ax = subplot(1, 1, 1);
            X = (1: size(image_orientation, 2)); % pixels
            Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_orientation);
            shading flat;
            colormap hsv;
            set(gca, 'YDir', 'reverse');
            caxis([0 180]);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
        end
        
        function obj = extract_local_orientation_3D(obj, PropertyName, wavelength, orientation, ratios, max_layer, SFB, SAR, K)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ratios: the depth in each ply
            % max_layer: the total plies of the structure
            % SFB, SAR: SpatialFrequencyBandwidth and SpatialAspectRatio
            % K: the smoothing applied to the Gabor magnitude responses. 
            % ***
            gaborArray             = gabor(wavelength, orientation, ...
                'SpatialFrequencyBandwidth', SFB, 'SpatialAspectRatio', SAR);
            wl_max                 = max(wavelength);
            Inplane_direction      = NaN(size(obj.(PropertyName)));
            Inp_direc_domain       = NaN(size(obj.(PropertyName)));
            Inp_direc_domain_layer = NaN(size(obj.(PropertyName)));
            ratio_spacing          = (ratios(end) - ratios(1)) / (length(ratios) - 1);
            for i = 1:max_layer
                ID_sum_layer = zeros(length(wavelength), length(orientation));
                for j = 1:length(ratios)
                    ratio                                          = ratios(j);
                    [Cscan_inam, profile_layer, profile_layer_pre] = obj.define_Cscan_depthprofile_knn(i, ratio, PropertyName, max_layer);
%                     % fillna and smooth the profile
%                     profile_layer                                  = fx_inpaint_nans(profile_layer, 5);
%                     profile_layer_pre                              = fx_inpaint_nans(profile_layer_pre, 5);
                    profile_layer                                  = medfilt2(profile_layer, [7, 7]);
                    profile_layer_pre                              = medfilt2(profile_layer_pre, [7, 7]);
                    % Apply filters to C-scan image.
                    [gaborMagnitude, ~]                            = imgaborfilt(Cscan_inam, gaborArray);
                    % normalized the magnitude.
                    for k = 1:length(gaborArray)
                        BW                       = gaborArray(k).SpatialFrequencyBandwidth;
                        sigmaX                   = gaborArray(k).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                        sigmaY                   = sigmaX ./ gaborArray(k).SpatialAspectRatio;
                        gabormag_k               = gaborMagnitude(:, :, k)  / (2 * sigmaX * sigmaY * pi);
                        % Gaussian filter for the GaborMag
                        sigma                    = 0.5*gaborArray(k).Wavelength;
                        invalid_part             = isnan(gabormag_k);
                        gabormag_k(invalid_part) = mean(gabormag_k, 'all', 'omitnan');
                        gabormag_k               = imgaussfilt(gabormag_k, K*sigma);
                        gabormag_k(invalid_part) = NaN;
                        gaborMagnitude(:, :, k)  = gabormag_k;
                    end
                    % extract the local orientaion
                    [lx, ly] = size(gaborMagnitude(:, :, 1));
                    %                     orientation_ave   = 0;
                    %                     count             = 0;
                    ID_sum            = zeros(length(wavelength), length(orientation));      
                    for ii = round(wl_max / 2): lx - round(wl_max)/2
                        for jj = round(wl_max / 2): ly - round(wl_max)/2
                            ID               = squeeze(gaborMagnitude(ii, jj, :));
                            if isnan(ID(1, 1)) % invalid part has been assigned NaN
                                Inplane_direction(ii, jj, k1:k2) = NaN;
                                continue;
                            end
                            ID          = reshape(ID, [obj.frow, obj.fcol]);
                            %                     ID_shift = circshift(ID, round(size(ID, 2) / 2), 2);
                            %                     ID_sub = ID - ID_shift;
                            [~, I]      = max(ID(:));
                            [~, indAng] = ind2sub(size(ID),I);
                            % sum up ID
                            ID_sum      = ID_sum + ID;
                            % save the orientation image
                            if ~isnan(profile_layer_pre(ii, jj)) && ...
                                    ~isnan(profile_layer(ii, jj)) && ...
                                    profile_layer_pre(ii, jj) > 0 && ...
                                    profile_layer_pre(ii, jj) <= profile_layer(ii, jj)
                                h                               = profile_layer(ii, jj) - profile_layer_pre(ii, jj);
                                k1                              = round(profile_layer_pre(ii, jj) + h*(ratio - ratio_spacing));
                                k2                              = round(profile_layer_pre(ii, jj) + h*(ratio + ratio_spacing));
                                Inp_direc_domain(ii, jj, k1:k2) = orientation(indAng);
                                %                                 % count the dominating angle (average)
                                %                                 orientation_ave                  = orientation_ave + orientation(indAng);
                                %                                 count                            = count + 1;
                            end
                        end
                    end
                    % analyze ID_sum
                    [~, I]                          = max(ID_sum(:));
                    [~, indAng]                     = ind2sub(size(ID_sum),I);
                    % assign the overall orientation
                    for ii = round(wl_max / 2): lx - round(wl_max)/2
                        for jj = round(wl_max / 2): ly - round(wl_max)/2
                            if ~isnan(profile_layer_pre(ii, jj)) && ...
                                    ~isnan(profile_layer(ii, jj)) && ...
                                    profile_layer_pre(ii, jj) > 0 && ...
                                    profile_layer_pre(ii, jj) <= profile_layer(ii, jj)
                                h                                = profile_layer(ii, jj) - profile_layer_pre(ii, jj);
                                k1                               = round(profile_layer_pre(ii, jj) + h*(ratio - ratio_spacing));
                                k2                               = round(profile_layer_pre(ii, jj) + h*(ratio + ratio_spacing));
                                Inplane_direction(ii, jj, k1:k2) = orientation(indAng);
                                % Inp_direc_domain(ii, jj, k1:k2) = orientation_ave / count;
                                % % count the dominating angle (average)
                                % orientation_ave                  = orientation_ave + orientation(indAng);
                                % count                            = count + 1;
                            end
                        end
                    end
                    ID_sum_layer      = ID_sum_layer + ID_sum;
                end
                % analyze ID_sum_layer
                [~, I]                          = max(ID_sum(:));
                [~, indAng]                     = ind2sub(size(ID_sum),I);
                % assign the layer orientation in one layer
                for ii = round(wl_max / 2): lx - round(wl_max)/2
                    for jj = round(wl_max / 2): ly - round(wl_max)/2
                        if ~isnan(profile_layer_pre(ii, jj)) && ...
                                ~isnan(profile_layer(ii, jj)) && ...
                                profile_layer_pre(ii, jj) > 0 && ...
                                profile_layer_pre(ii, jj) <= profile_layer(ii, jj)
                            h                                = profile_layer(ii, jj) - profile_layer_pre(ii, jj);
                            k1                               = round(profile_layer_pre(ii, jj) + h*(ratio - ratio_spacing));
                            k2                               = round(profile_layer_pre(ii, jj) + h*(ratio + ratio_spacing));
                            Inp_direc_domain_layer(ii, jj, profile_layer_pre(ii, jj):profile_layer(ii, jj)) = orientation(indAng);
                            % Inp_direc_domain(ii, jj, k1:k2) = orientation_ave / count;
                            % % count the dominating angle (average)
                            % orientation_ave                  = orientation_ave + orientation(indAng);
                            % count                            = count + 1;
                        end
                    end
                end
                % end loop in one Cscan
                disp([num2str(i), '/', num2str(max_layer), ' layer']);
            end
            obj.Inplane_direction_3D               = Inplane_direction;
            obj.Inplane_direction_3D_overall       = Inp_direc_domain;
            obj.Inplane_direction_3D_overall_layer = Inp_direc_domain_layer;
        end
        
        function obj = extract_local_orientation_3D_parallel_allwl(obj, PropertyName, wavelength, orientation, ratios, K, sigma_denoise)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ratios: the ratio determing the index by the front and back surfaces , unit: arb.
            % K: the smoothing applied to the Gabor magnitude responses.
            % sigma: threshould for curvelet denoising, 0 means no application of denoising
            % ***
            obj.wavelength_LG         = wavelength;
            obj.theta_LG              = orientation;
            obj.PropertyName_LG       = PropertyName;       
            % start 
            gaborArray                = gabor(wavelength, orientation, ...
                'SpatialFrequencyBandwidth', 1, 'SpatialAspectRatio', 0.7);
            row_ga                    = length(wavelength);
            col_ga                    = length(orientation);
            wl_max                    = round(max(wavelength)/2);
            index_container           = NaN(size(obj.front_I));
            Inplane_direction         = NaN(size(obj.(PropertyName)));
            Inplane_direction_overall = NaN(size(obj.(PropertyName)));
            ID_sum_all = NaN(length(ratios), length(wavelength), length(orientation)); 
            %
            disp('gaborMagnitude calculation...');
            for i = 1:length(ratios)
                ratio = ratios(i); % the ratio to determine the index
                if i~=length(ratios)
                    ratio_next = ratios(i+1);
                end
                index_container_prev = index_container;
                [~, C_scan_inam_para, index_container] = obj.define_parallel_inamCscan(ratio, PropertyName);
                [~, ~, index_container_next]           = obj.define_parallel_inamCscan(ratio_next, PropertyName);
                if sigma_denoise==0
                    C_scan_inam_para_denoise = C_scan_inam_para;
                else
                    C_scan_inam_para_denoise = fx_curvelet_denoise_enhanced(sigma_denoise, C_scan_inam_para);
                end
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
                        if i==1
                            upper_index_bound = index_container(ii, jj);
                            lower_index_bound = round(index_container(ii, jj)/2 ...
                            +index_container_next(ii, jj)/2);
                        elseif i==length(ratios)
                            upper_index_bound = round(index_container(ii, jj)/2 ...
                                + index_container_prev(ii, jj)/2);
                            lower_index_bound = index_container(ii, jj);
                        else
                            upper_index_bound = round(index_container(ii, jj)/2 ...
                                + index_container_prev(ii, jj)/2);
                            lower_index_bound = round(index_container(ii, jj)/2 ...
                                +index_container_next(ii, jj)/2);
                        end
                        Inplane_direction(ii, jj, upper_index_bound: ...
                            lower_index_bound) = orientation(indAng);
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
                disp([num2str(ratios(i)), '/', num2str(ratios(end))]);
            end
            obj.Inplane_direction_3D_ID         = Inplane_direction;
            obj.Inplane_direction_3D_overall_ID = Inplane_direction_overall;
            obj.ID_sum_overall_ID               = ID_sum_all;
        end
        
        function obj = extract_local_orientation_3D_plywise_allwl(obj, PropertyName, wavelength, orientation, ratios, plies, K, sigma_curvelet)
            % extract the local orientation image in each ply
            % adjust the scale by the edge
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ratios: the ratios to define the plane in one ply
            % plies: number of the plies
            % K: the smoothing applied to the Gabor magnitude responses.
            % sigma_curvelet: threshould for curvelet denoising, 0 means no application of denoising
            % ***
            obj.wavelength_LG         = wavelength;
            obj.theta_LG              = orientation;
            obj.PropertyName_LG       = PropertyName;
            % start
            gaborArray                = gabor(wavelength, orientation, ...
                'SpatialFrequencyBandwidth', 1, 'SpatialAspectRatio', 0.5);
            row_ga                    = length(wavelength);
            col_ga                    = length(orientation);
            wl_max                    = max(wavelength);
            wl_min                    = min(wavelength);
            index_container           = NaN(size(obj.front_I));
            Inplane_direction         = NaN(size(obj.(PropertyName)));
            Inplane_direction_overall = NaN(size(obj.(PropertyName)));
            ID_sum_all                = NaN(plies, length(ratios), ...
                length(wavelength), length(orientation));
            [lx, ly]                  = size(obj.front_I);
            disp('gaborMagnitude calculation...');
            for p = 1:plies
                tic;
                for ratio_index = 1:length(ratios)
                    %                     tic;
                    ratio                                 = ratios(ratio_index);
                    index_container_prev                  = index_container;
                    [~, C_scan_inam_ply, index_container] = obj.define_plywise_inamCscan(p, ratio, PropertyName);
                    if ratio_index==1
                        [~, ~, upper_index_bound]    = obj.define_plywise_inamCscan ...
                            (p, 0, PropertyName);
                        [~, ~, index_container_next] = obj.define_plywise_inamCscan ...
                            (p, ratios(ratio_index+1), PropertyName);
                        lower_index_bound            = round(index_container/2 ...
                            +index_container_next/2);
                    elseif ratio_index==length(ratios)
                        upper_index_bound            = round(index_container/2 ...
                            + index_container_prev/2);
                        [~, ~, lower_index_bound]    = obj.define_plywise_inamCscan ...
                            (p, 1, PropertyName);
                    else
                        [~, ~, index_container_next] = obj.define_plywise_inamCscan ...
                            (p, ratios(ratio_index+1), PropertyName);
                        upper_index_bound = round(index_container/2 ...
                            + index_container_prev/2);
                        lower_index_bound = round(index_container/2 ...
                            +index_container_next/2);
                    end  
                    % 2D log-gabor filter + ID multi-resolution analysis
                    if sigma_curvelet==0
                        % fill nan
                        C_scan_inam_para_denoise = fx_inpaint_nans(C_scan_inam_ply, 5);
%                         C_scan_inam_para_denoise = C_scan_inam_ply;
                    else
                        C_scan_inam_para_denoise = fx_curvelet_denoise_enhanced(sigma_curvelet, C_scan_inam_ply);
                    end
                    [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_para_denoise, gaborArray);
                    % normalized the magnitude.
                    for k = 1:length(gaborArray)
                        BW                      = gaborArray(k).SpatialFrequencyBandwidth;
                        sigmaX                  = gaborArray(k).Wavelength / pi * sqrt(log(2)/2) * (2^BW+1) / (2^BW-1);
                        sigmaY                  = sigmaX ./ gaborArray(k).SpatialAspectRatio;
                        % normalized the magnitude.
                        gabormag_k              = gaborMagnitude(:, :, k);
                        gabormag_k              = gabormag_k / (2 * sigmaX * sigmaY * pi);
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
%                     toc;
                    % search for the max
                    ID_sum         = zeros(length(wavelength), length(orientation));
                    % check the nan and adjust the ID
                    NaN_part       = isnan(C_scan_inam_ply);
                    [r_nan, l_nan] = find(NaN_part);
                    % calculate the matrix of the distances to nan part
                    [Y, X]         = meshgrid(1:ly, 1:lx);
                    distances_pow_2 = max(wavelength)^2 * ones(lx, ly);
                    for i = 1:length(r_nan)
                        distances_i_pow2 = (r_nan(i)-X).^2 + (l_nan(i)-Y).^2;
                        distances_pow_2   = min(distances_pow_2, distances_i_pow2);
                    end
                    for ii = round(wl_min) + 1: lx - round(wl_min)-1
                        for jj = round(wl_min) + 1: ly - round(wl_min)-1
                            if obj.rear_I(ii, jj) < index_container(ii, jj)% check rear surface
                                continue;
                            end
                            ID = squeeze(gaborMagnitude(ii, jj, :));
                            ID = reshape(ID, [row_ga, col_ga]);
                            % check the nan and adjust the ID
                            %                             distances = sqrt((r_nan-ii).^2 + (l_nan-jj).^2);
                            % binary Search the first element that are
                            % larger than the distances index from 1 to max+1
                            %                             ID_bound  = fx_binarySearch(...
                            %                                 wavelength, length(wavelength), min(distances));
                            % wavelength is min:2:max!!
                            %                             ID_bound = min(round(min(distances)), wavelength(end));
                            ID_bound = round(sqrt(distances_pow_2(ii, jj)));
                            if ID_bound <= 4
                                continue;
                            elseif ID_bound < wl_max || ... % consider the edge of the whole image
                                    ii < wl_max || jj < wl_max  || lx-ii < wl_max || ly-jj < wl_max 
                                % consider the edge of the whole image
                                ID_bound     = min([ID_bound, ii, jj, lx-ii , ly-jj]);
                                ID_bound_ind = fx_binarySearch(wavelength, length(wavelength), ID_bound)-1;
                                ID           = ID(1:ID_bound_ind, :);
                            end
                            [~, I]      = max(ID(:));
                            [~, indAng] = ind2sub(size(ID),I);
                            if ~isnan(upper_index_bound(ii, jj)) && ...
                                    ~isnan(lower_index_bound(ii, jj))
                                Inplane_direction(ii, jj, upper_index_bound(ii, jj): ...
                                    lower_index_bound(ii, jj)) = orientation(indAng);
                            end
                            % sum up the ID
                            if size(ID)==size(ID_sum)
                                ID_sum = ID_sum + ID / max(ID, [], 'all');
                            end
                        end
                    end
                    % save the summary ID
                    ID_sum_all(p, ratio_index, :, :) = ID_sum;
                    %                     [~, I]                                 = max(ID_sum(:));
                    %                     [~, indAng]                            = ind2sub(size(ID_sum),I);
                    %                     Inplane_direction_overall(:, :, upper_index_bound: ...
                    %                         lower_index_bound) = orientation(indAng);
%                     toc;
%                     disp([num2str(p), ' / ', num2str(plies), ':', num2str(ratio_index)]);
                end
                %                 % end loop in one Cscan
                toc;
                disp([num2str(p), ' / ', num2str(plies)]);
            end
            obj.Inplane_direction_3D_ID         = Inplane_direction;
            obj.Inplane_direction_3D_overall_ID = Inplane_direction_overall;
            obj.ID_sum_overall_ID               = ID_sum_all;
        end
        
        function obj = extract_local_ID_AlongDepth_LG(obj, PropertyName, wavelength, orientation, distances_to_front, folder)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength, orientation: the ratio for the multi-resolution Gabor filter bank
            % folder: the foldname of storage
            % ***
            gaborArray = gabor(wavelength, orientation);
            inam = abs(obj.(PropertyName));
            [lx, ly, ~] = size(inam);
            % extract the local orientaion
            lgrow = length(wavelength);
            lgcol = length(orientation);
            for d = 1:length(distances_to_front)
                dis_to_front = distances_to_front(d) * obj.fs / 1e6; % us to sample points
                C_scan_inam_para = NaN([lx, ly]);
                for ii = 1: size(inam, 1)
                    for jj = 1:size(inam, 2)
                        front_index = obj.front_I(ii, jj);
                        if front_index + dis_to_front < obj.rear_I(ii, jj) && front_index > 0
                            C_scan_inam_para(ii, jj) = inam(ii, jj, round(front_index + dis_to_front));
                        end
                    end
                end
                % Radon transform
                % fillna
                C_scan_inam_para(isnan(C_scan_inam_para)) = mean(C_scan_inam_para, 'all', 'omitnan');
                % pad the C_scan_inam matrix
                wl_max = max(wavelength) / 2;
                C_scan_inam_RT_pad = padarray(C_scan_inam_para, [round(wl_max) round(wl_max)], ...
                    mean(C_scan_inam_para, 'all', 'omitnan'));
                [gaborMag_d, ~] = imgaborfilt(C_scan_inam_RT_pad, gaborArray);
                % normalized the magnitude.
                for i = 1:length(gaborArray)
                    BW = gaborArray(i).SpatialFrequencyBandwidth;
                    sigmaX = gaborArray(i).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                    sigmaY = sigmaX ./ gaborArray(i).SpatialAspectRatio;
                    gaborMag_d(:, :, i) = gaborMag_d(:, :, i) / (2 * sigmaX * sigmaY * pi);
                end
                % create all the files and write basic info.
                jj = round(ly / 2);
                if d==1
                    mkdir(folder);
                    for ii = 1: lx
                        baseFileName = ['x_' num2str(ii), '_', 'y_', num2str(jj), '_ID.txt'];
                        fullFileName = fullfile(folder, baseFileName);
                        [fileID, ~]  = fopen(fullFileName, 'w');
                        fprintf(fileID,'%6s \n', 'lgrow');
                        fprintf(fileID,'%6f \n', lgrow);
                        fprintf(fileID,'%6s \n', 'lgcol');
                        fprintf(fileID,'%6f \n', lgcol);
                        fclose(fileID);
                    end
                end
                %
                for ii = 1: lx
                    if obj.front_I(ii, jj) + dis_to_front < obj.rear_I(ii, jj)
                        ID_d         = squeeze(gaborMag_d(ii + round(wl_max), jj + round(wl_max), :));
                        ID_shape     = reshape(ID_d, [lgrow, lgcol]);
                        baseFileName = ['x_' num2str(ii), '_', 'y_', num2str(jj), '_ID.txt'];
                        fullFileName = fullfile(folder, baseFileName);
                        [fileID, ~]  = fopen(fullFileName, 'a');
                        fprintf(fileID,'%6s \n', ' ');
                        fprintf(fileID,'%6.2f \n', dis_to_front);
                        fprintf(fileID,'%6f \n', ID_shape);
                        fclose(fileID);
                    end
                end
                disp([num2str(distances_to_front(d)) '/' num2str(distances_to_front(end))]);
            end
            %             obj.ID_pixel = ID;
            %             % display the ID along the depth for each r_idx
            %             cf = figure('Name', ['ID' '_' 'LG' '_' 'x' num2str(x) '_' 'y' num2str(y)]);
            %             set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            %             ax = subplot(1, 1, 1);
            %             imagesc(ax, angular_1D);
            %             colormap jet;
            %             h = colorbar;
            %             set(get(h, 'Title'), 'string', ...
            %                 '\fontname {times new roman}\fontsize {16} Mag. of angular distribution from RT (arb.)');
            %             set(ax, 'fontsize', 16);
            %             set(ax, 'linewidth', 1.5);
            %             xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Depth (pixel)', 'fontsize', 16);
            %             set(gca, 'ZDir', 'reverse');
        end
        
        function show_inplane_direction_3D_ID(obj, xslice, yslice, zslice, mfsize, angle_compens)
            % demonstrate the 3d inplane_direction extracted by Gabor filter
            % demonstrate the ID_sum_overall
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % the properties to show are: obj.Inplane_direction_3D_ID,
            % obj.Inplane_direction_3D_overall_ID,
            % obj.ID_sum_overall_ID;
            % angle_compens: compensation for angle
            % mfsize: the medfilter size
            % ***
            inph_ex = obj.Inplane_direction_3D_ID;
            inph_ex = mod(inph_ex + angle_compens, 180) - 90;
            inph_ex = medfilt3(inph_ex, mfsize);
            y       = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x       = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            % ******************** show ply tracks
            cf = figure('Name', ['3d_orientation_2Dfilter_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            x_idx = xslice * obj.fx / 1e3;
            y_idx = yslice * obj.fy / 1e3;
            mask  = (obj.col==x_idx) | (obj.row==y_idx);
            x_row = obj.row(mask) / obj.fx * 1e3;
            x_col = obj.col(mask) / obj.fy * 1e3;  % / obj.fy * 1e3 ;
            x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
            if sum(isnan(obj.front_I))~=0
                h1 = scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                    y, obj.front_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
                    3, 'black', 'filled', 'DisplayName','Front surface');
                hold on;
                h2 = scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                    y, obj.rear_I(:, x_idx)/ obj.fs * 1e3 * 3000/2, ...
                    3, [0.75 0.75 0.75], 'filled', 'DisplayName','Rear surface');
                hold on;
                h3 = scatter3(ax, x_col, x_row, x_dep, ...
                    3, [122 122 121]/255, 'filled', ...
                    'DisplayName','Interply track');
                hold on;
                % select the yslice
                scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                    obj.front_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
                    3, 'black', 'HandleVisibility','off');
                hold on;
                scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                    obj.rear_I(y_idx, :)/ obj.fs * 1e3 * 3000/2, ...
                    3, [0.75 0.75 0.75], 'HandleVisibility','off');
                hold on;
            end
            % the 1st dimension is length(radius)
            if size(size(inph_ex))==4  % check dims
                z         = (0: size(inph_ex, 4) - 1) / obj.fs * 1e3 * 3000/2;
                [X, Y, Z] = meshgrid(x, y, z);
                % inph_visual(mask_phase_interp_visual==1) = NaN;
                for i = 1: size(inph_ex, 1)
                    h  = slice(ax, X, Y, Z, squeeze(inph_ex(i,:, :, :)), xslice, yslice, zslice);
                    hold on;
                    set(h, 'EdgeColor', 'none');
                    colormap hsv;
                    h  = colorbar;
                    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
                    set(ax, 'fontsize', 16);
                    set(ax, 'linewidth', 1.5);
                    xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
                    ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
                    zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
                    set(gca, 'ZDir', 'reverse');
                end
            else
%                 z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e3 * 3000/2;
                z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6; % time axis, unit: us
                [X, Y, Z] = meshgrid(x, y, z);
                h         = slice(ax, X, Y, Z, inph_ex , xslice, yslice, zslice);
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
            end
            %             legend([h1 h2 h3], 'Front surface', 'Rear surface', 'Interply track');
            xlim([x(1) x(end)]);
            ylim([y(1) y(end)]);
            %             legend([h3 h1 h2], 'Interply track', 'Front surface', 'Rear surface');
            %             ** Uncommnet for non-interply demp!
            legend([h1 h2], 'Front surface', 'Rear surface');
            xlim([0 x(end)]);
            ylim([0 y(end)]);
%             zlim([0.8 7]);
%             zticks([0 1 2 3 4 5 6]);
            view([15 65 40]);
%             zlim([0 3]);
%             zlim([1 5.5])
%             view([15 45 30]);
%              % ************ Inplane_orientation_3D_overall ******
%             inph_ex   = obj.Inplane_direction_3D_overall_ID;
%             y         = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
%             x         = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
%             z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
%             [X, Y, Z] = meshgrid(x, y, z);
%             cf        = figure('Name', ['3d_orientation_RT_overall_' num2str(obj.radius_RT) '_' num2str(xslice(1))]);
%             set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
%             ax        = subplot(1, 1, 1);
%             h         = slice(ax, X, Y, Z, inph_ex, xslice, yslice, zslice);
%             hold on;
%             set(h, 'EdgeColor', 'none');
%             colormap hsv;
%             h = colorbar;
%             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (degree)');
%             set(ax, 'fontsize', 16);
%             set(ax, 'linewidth', 1.5);
%             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
%             ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
%             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
%             set(gca, 'ZDir', 'reverse');
%             % ************ Inplane_orientation_3D_overall ******
%             ID_sum_all     = obj.ID_sum_overall_ID;
%             for i = 1:size(ID_sum_all, 1)
%                 ID_sum_all(i, :) = ID_sum_all(i, :) / max(ID_sum_all(i, :));
%             end
%             cf = figure('Name', ['3d_ID_sum_overall' '_' num2str(obj.radius_RT)]);
%             set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
%             ax = subplot(1, 1, 1);
%             imagesc(ax, linspace(0,180,size(ID_sum_all, 2)), 1:size(ID_sum_all, 1), squeeze(sum(ID_sum_all, 2)));
%             hold on;
%             set(h, 'EdgeColor', 'none');
%             colormap jet;
%             h = colorbar;
%             set(get(h, 'Title'), 'string', ...
%                 '\fontname {times new roman}\fontsize {16} Mag. of angular distribution from RT (arb.)');
%             set(ax, 'fontsize', 16);
%             set(ax, 'linewidth', 1.5);
%             xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
%             ylabel('\fontname {times new roman} Depth (pixel)', 'fontsize', 16);
%             set(gca, 'ZDir', 'reverse');
        end
        
        function [m_fiber_angle, std_fiber_angle, yhat, Idof_N] = calculate_m_std_fiber_angle(obj, PropertyName, PropertyName_IFD, ratio, ratio_next, angle_compens)
            % calculate the mean fiber angle and its standard deviation
            % PropertyName: the property for defining the index
            % PropertyName_IFD: the property for extracting the inplane fiber direction
            % ratio, ratio_next: the ratios for defining the index aligned with the surfaces
            % ref_angle: the ideal fiber angle in the ply, used for calculation
            [~, ~, index_container]  = obj.define_parallel_inamCscan(ratio, PropertyName);
%             [~, ~, index_container_next] = obj.define_parallel_inamCscan(ratio_next, PropertyName);
            [lx, ly] = size(index_container);
            Inplane_direction        = obj.(PropertyName_IFD) + angle_compens;
%             Inplane_direction_oneply     = Inplane_direction;
            Inplane_direction_oneply = NaN(lx, ly);
            for ii = 1: lx
                for jj = 1: ly
                    % upper_index_bound = round(index_container(ii, jj));
                    % lower_index_bound = round(index_container_next(ii, jj));
                    % Inplane_direction_oneply(ii, jj, 1:upper_index_bound) = NaN;
                    % Inplane_direction_oneply(ii, jj, lower_index_bound:end) = NaN;
                    Inplane_direction_oneply(ii, jj) = Inplane_direction(...
                        ii, jj, round(index_container(ii, jj)));
                end
            end
            % minus the ref_angle to get the diff
            %             Inplane_direction_oneply = Inplane_direction_oneply - ref_angle;
            Idof_edge    = [0:180] - 22;
            Idof_N       = histcounts(Inplane_direction_oneply, Idof_edge);
%             Idof_N       = circshift(Idof_N, 22); % shift 22 degree for fitting
            %************************* Gaussian fitting ********************
            centers      = [0 45 90 135] -22;
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
            std_fiber_angle = amp(amp_I);
            % *************************
            % make sure the diff is less than 90 or larger than -90
            %             Inplane_direction_oneply(Inplane_direction_oneply>90) = ...
            %                 Inplane_direction_oneply(Inplane_direction_oneply>90)- 180;
            %             Inplane_direction_oneply(Inplane_direction_oneply<-90) = ...
            %                 Inplane_direction_oneply(Inplane_direction_oneply<-90)+ 180;
            %             % set threshold for cu
            %             Inplane_direction_oneply(Inplane_direction_oneply>45)  = NaN;
            %             Inplane_direction_oneply(Inplane_direction_oneply<-45) = NaN;
            %
%             std_fiber_angle = parameter(amp_I*2);
%             m_fiber_angle_diff   = mean(Inplane_direction_oneply, 'all','omitnan');
%             m_fiber_angle        = m_fiber_angle_diff + ref_angle;
            std_fiber_angle      = std(Inplane_direction_oneply, 0, 'all','omitnan');
        end
        
        function Idof_N = calculate_angle_distrib(obj, PropertyName, PropertyName_IFD, ratio, angle_compens)
            % calculate the mean fiber angle and its standard deviation
            % PropertyName: the property for defining the index
            % PropertyName_IFD: the property for extracting the inplane fiber direction
            % ratio, ratio_next: the ratios for defining the index aligned with the surfaces
            % ref_angle: the ideal fiber angle in the ply, used for calculation
            [~, ~, index_container]  = obj.define_parallel_inamCscan(ratio, PropertyName);
            %             [~, ~, index_container_next] = obj.define_parallel_inamCscan(ratio_next, PropertyName);
            [lx, ly] = size(index_container);
            Inplane_direction        = obj.(PropertyName_IFD);
            %             Inplane_direction_oneply     = Inplane_direction;
            Inplane_direction_oneply = NaN(lx, ly);
            for ii = 1: lx
                for jj = 1: ly
                    % upper_index_bound = round(index_container(ii, jj));
                    % lower_index_bound = round(index_container_next(ii, jj));
                    % Inplane_direction_oneply(ii, jj, 1:upper_index_bound) = NaN;
                    % Inplane_direction_oneply(ii, jj, lower_index_bound:end) = NaN;
                    Inplane_direction_oneply(ii, jj) = Inplane_direction(...
                        ii, jj, round(index_container(ii, jj)));
                end
            end
            % minus the ref_angle to get the diff
            %             Inplane_direction_oneply = Inplane_direction_oneply - ref_angle;
            Idof_edge = [0.5:180.5];
            Idof_N    = histcounts(Inplane_direction_oneply, Idof_edge);
            Idof_N    = circshift(Idof_N, 22+angle_compens);
        end
        
        function show_inplane_direction_3D_ID_Bscan(obj, B_type, index, win, mfsize, angle_compens)
            % demonstrate the B_scan of the interply track
            % show inph_ex
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % win: the win of the B-scan image
            % angle_compens: compensation angle
            % mfsize: the medfilter size
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            inph_ex = obj.Inplane_direction_3D_ID;
            inph_ex = mod(inph_ex + angle_compens, 180);
            inph_ex = medfilt3(inph_ex, mfsize);
            x           = (1:(win(end)-win(1)+1))/obj.fx*1e3;
            z           = (1:size(inph_ex, 3))/obj.fs*1e6;
            cf          = figure('Name', ['B_scan_Inplane_angles_', B_type, '_', num2str(index)]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            if (B_type == 'x')
                % B scan image
                B_scan      = real(squeeze(inph_ex(index, win, :)));
                %                 imagesc(x, z, B_scan');
                pcolor(x, z, B_scan');
                shading flat;
                set(gca, 'ydir', 'reverse');
                hold on;
                mask  = (obj.row==index);
                x_row = obj.col(mask) / obj.fx * 1e3;
                x_dep = obj.dep(mask) / obj.fs * 1e6;  % / obj.fy * 1e3 ;
                h1 = scatter(ax, x_row, x_dep, ...
                    3,  [122 122 121]/255, 'filled', ...
                    'DisplayName','Interply track');
                hold on;
                h2 = scatter(x, obj.front_I(index, :)/ obj.fs * 1e6, ...
                    3, 'red', 'filled', ...
                    'DisplayName','Front surface');
                hold on;
                h3 = scatter(x, obj.rear_I(index, :)/ obj.fs * 1e6, ...
                    3, 'magenta', 'filled', ...
                    'DisplayName','Rear surface');
                hold on;
                % figure setup
                colormap(hsv);
                h = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle. (degree)');
                caxis([0 180]);
                %             cl   = caxis;
                %             caxis(cl);
                set(gca, 'Fontname', 'times new Roman');
                xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
                ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            elseif (B_type == 'y')
                % B scan image
                B_scan      = real(squeeze(inph_ex(win, index, :)));
                %                 imagesc(x, z, B_scan');
                pcolor(x, z, B_scan');
                shading flat;
                set(gca, 'ydir', 'reverse');
                hold on;
                mask  = (obj.col==index);
                x_row = obj.row(mask) / obj.fy * 1e3;
                x_dep = obj.dep(mask) / obj.fs * 1e6;  % / obj.fy * 1e3 ;
                h1 = scatter(ax, x_row, x_dep, ...
                    3,  [122 122 121]/255, 'filled', ...
                    'DisplayName','Interply track');
                hold on;
                h2 = scatter(x, obj.front_I(:, index)/ obj.fs * 1e6, ...
                    3, 'red', 'filled', ...
                    'DisplayName','Front surface');
                hold on;
                h3 = scatter(x, obj.rear_I(:, index)/ obj.fs * 1e6, ...
                    3, 'magenta', 'filled', ...
                    'DisplayName','Rear surface');
                hold on;
                % figure setup
                colormap(hsv);
                h = colorbar;
                set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle. (degree)');
                %             cl   = caxis;
                %             caxis(cl);
                xlabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
                ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            end
            legend([h2 h3 h1], 'Front surface', 'Rear surface', 'Interply track');
            set(ax, 'fontsize', 16, 'Fontname', 'times new Roman');
            set(ax, 'linewidth', 2);
        end
        
        function show_inplane_direction_3D_ID_meanstd(obj, xslice, yslice, zslice, plies, ratios, PropertyName, wavelength)
            % demonstrate the 3d inplane_direction extracted by Gabor filter
            % demonstrate the ID_sum_overall
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % ************ mean ************
            inph_ex                         = obj.Inplane_direction_3D_ID;
            [lx, ly, ~]                     = size(inph_ex);
            Inplane_d_plywise_mean          = NaN(size(inph_ex));
            Inplane_d_plywise_std           = NaN(size(inph_ex));
            Inplane_d_plywise_mean_wholeply = NaN(size(inph_ex));
            wl_max                          = max(wavelength);
            plywise_direction               = NaN(size(inph_ex, 1), size(inph_ex, 2), plies, length(ratios));
            for p = 1:plies
                tic;
                % check the mean and std values of the orientations in one ply
                [~, ~, upper_index_bound]    = obj.define_plywise_inamCscan ...
                    (p, 0, PropertyName);
                [~, ~, lower_index_bound]    = obj.define_plywise_inamCscan ...
                    (p, 1, PropertyName);
                for ratio_index = 1:length(ratios)
                    for ii = round(wl_max / 2) + 1: round(lx - wl_max / 2)
                        for jj = round(wl_max / 2) + 1: round(ly - wl_max / 2)
                            plywise_direction(ii, jj, p, ratio_index) = mode(inph_ex(ii, jj, upper_index_bound(ii, jj): ...
                                lower_index_bound(ii, jj)), 'all');
                            Inplane_d_plywise_mean(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj)) ...
                                = mode(inph_ex(ii, jj, upper_index_bound(ii, jj): ...
                                lower_index_bound(ii, jj)), 'all');
                            Inplane_d_plywise_std(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj)) ...
                                = sum( abs(inph_ex(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj)) - ...
                                Inplane_d_plywise_mean(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj))), 'all', 'omitnan');
                        end
                    end
                end
                % assign the orientation for a ply
                for ii = round(wl_max / 2) + 1: round(lx - wl_max / 2)
                    for jj = round(wl_max / 2) + 1: round(ly - wl_max / 2)
                        Inplane_d_plywise_mean_wholeply(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj)) ...
                            = mode(plywise_direction(:, :, p, :), 'all');
                    end
                end
                toc;
            end
            % ********* display mean **************
            y         = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x         = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            cf        = figure('Name', ['3d_orientation_2Dfilter_mean_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax        = subplot(1, 1, 1);
            h         = slice(ax, X, Y, Z, Inplane_d_plywise_mean, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap hsv;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mean angle (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            % ********* std **************
            y         = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x         = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            cf        = figure('Name', ['3d_orientation_2Dfilter_mean_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax        = subplot(1, 1, 1);
            h         = slice(ax, X, Y, Z, Inplane_d_plywise_std , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} STD (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            % ********* mean whole ply **************
            y         = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x         = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z         = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            cf        = figure('Name', ['3d_orientation_2Dfilter_mean_wholeply_', 'xslice', num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax        = subplot(1, 1, 1);
            h         = slice(ax, X, Y, Z, Inplane_d_plywise_mean_wholeply , xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h         = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} STD (degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
        end
        
        function show_inplane_direction_3D_ID_global(obj, plies, ratios, orientation)
            % demonstrate the global inplane_direction by the 'ply ID'
            ID_sum_all               = obj.ID_sum_overall_ID;              
            L_ra                     = length(ratios);       
            Inplane_direction_global = NaN(1, L_ra);
            for p = 1:plies
                for ra_index = 1:L_ra
                    ID                                            = squeeze(ID_sum_all(p, ra_index, :, :));
                    [~, I]                                        = max(ID(:));
                    [~, indAng]                                   = ind2sub(size(ID),I);
                    Inplane_direction_global(L_ra*(p-1)+ra_index) = orientation(indAng);
                end
            end
            figure,
            scatter(1:length(Inplane_direction_global), ...
                Inplane_direction_global, 10, 'filled');
            hold on;
            for p = 1:plies
                xline(p * L_ra);
            end
        end
        
        function save_local_orientaion_3D(obj, savename)
            % save all the file of the 3D orientation
            % savename: the name incl. the path and the file for saving
            % save member variables
            Inplane_direction_3D_var = obj.Inplane_direction_3D;
            Inplane_direction_3D_overall_var = obj.Inplane_direction_3D_overall;
            ID_sum_overall_var = obj.ID_sum_overall;
            Inplane_direction_3D_ID_var = obj.Inplane_direction_3D_ID;
            Inplane_direction_3D_overall_ID_var = obj.Inplane_direction_3D_overall_ID;
            ID_sum_overall_ID_var = obj.ID_sum_overall_ID;
            Inplane_direction_3D_modifyRT_var = obj.Inplane_direction_3D_modifyRT;
            Inplane_direction_3D_multireso_modifyRT_var = obj.Inplane_direction_3D_multireso_modifyRT;
            Inplane_direction_3D_overall_modifyRT_var = obj.Inplane_direction_3D_overall_modifyRT;
            ID_sum_overall_modifyRT_var = obj.ID_sum_overall_modifyRT;
            % supplement info.
            radius_RT_var = obj.radius_RT;
            theta_RT_var = obj.theta_RT;
            distances_to_front_RT_var = obj.distances_to_front_RT;
            PropertyName_RT_var = obj.PropertyName_RT;
            wavelength_LG_var = obj.wavelength_LG;
            theta_LG_var = obj.theta_LG;
            distances_to_front_LG_var = obj.distances_to_front_LG;
            PropertyName_LG_var = obj.PropertyName_LG;
            save(savename,'Inplane_direction_3D_var',...
                'Inplane_direction_3D_overall_var', ...
                'ID_sum_overall_var', ...
                'Inplane_direction_3D_ID_var',...
                'Inplane_direction_3D_overall_ID_var', ...
                'ID_sum_overall_ID_var', ...
                'radius_RT_var', ...
                'theta_RT_var', ...
                'distances_to_front_RT_var', ...
                'PropertyName_RT_var', ...
                'wavelength_LG_var', ...
                'theta_LG_var', ...
                'distances_to_front_LG_var', ...
                'PropertyName_LG_var', ...
                'Inplane_direction_3D_modifyRT_var', ...
                'Inplane_direction_3D_multireso_modifyRT_var', ...
                'Inplane_direction_3D_overall_modifyRT_var', ...
                'ID_sum_overall_modifyRT_var',...
                '-v7.3');
        end
        
        function obj = load_local_orientaion_3D(obj, loadname)
            % save all the file of the 3D orientation
            % loadname: the name incl. the path and the file for loading
            % load member variables
            load(loadname);
            obj.Inplane_direction_3D                    = Inplane_direction_3D_var;
            obj.Inplane_direction_3D_overall            = Inplane_direction_3D_overall_var;
            obj.ID_sum_overall                          = ID_sum_overall_var;
            obj.Inplane_direction_3D_ID                 = Inplane_direction_3D_ID_var;
            obj.Inplane_direction_3D_overall_ID         = Inplane_direction_3D_overall_ID_var;
            obj.ID_sum_overall_ID                       = ID_sum_overall_ID_var;
            obj.Inplane_direction_3D_modifyRT           = Inplane_direction_3D_modifyRT_var;
            obj.Inplane_direction_3D_multireso_modifyRT = Inplane_direction_3D_multireso_modifyRT_var;
            obj.Inplane_direction_3D_overall_modifyRT   = Inplane_direction_3D_overall_modifyRT_var;
            obj.ID_sum_overall_modifyRT                 = ID_sum_overall_modifyRT_var;
            % supplement info.
            obj.radius_RT             = radius_RT_var;
            obj.theta_RT              = theta_RT_var;
            obj.distances_to_front_RT = distances_to_front_RT_var;
            obj.PropertyName_RT       = PropertyName_RT_var;
            obj.wavelength_LG         = wavelength_LG_var;
            obj.theta_LG              = theta_LG_var;
            obj.distances_to_front_LG = distances_to_front_LG_var;
            obj.PropertyName_LG       = PropertyName_LG_var;
        end
        
        function makeMovie_fiberDirection(obj, B_type, win, property_name, angle_compens)
            % make a movie of the B_scan of the interply track
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % win: the win of the B-scan image
            % property_name: name of the dataset
            % angle_compens: compensation angle
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            dataset = obj.(property_name);
            dataset = mod(dataset + angle_compens, 180) - 90;
            front_p = obj.front_I;
            bakc_p  = obj.rear_I;
            x       = (1:(win(end)-win(1)+1))/obj.fx * 1e3;
            y       = (1:(win(end)-win(1)+1))/obj.fy * 1e3;
            z       = (1:size(dataset, 3))/obj.fs * 1e3 * 3000/2;
            cf      = figure('Name', ['B_scan_angles_', B_type, '_', num2str(win(end))]);
            set(cf, 'Position', [0, 0, length(win)*2, 1200], 'color', 'white');
            % make vedio
            axis tight manual;
            title(['y =' num2str(0/obj.fy * 1e3) ' mm']);
            xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} z (mm)', 'fontsize', 16);
            v = VideoWriter([property_name, '.avi']);
            open(v);
            %
            if (B_type == 'x')
                for index = 1:size(dataset,1)
                    % B scan image
                    % side view
                    B_scan = squeeze(dataset(index, win, :));
                    figure(cf);
                    ca1 = subplot(2, 1, 1);
                    % imagesc(x, z, B_scan');
                    pcolor(x, z, B_scan');
                    shading flat;
                    set(ca1, 'ydir', 'reverse');
                    hold on;
                    title(['y =' num2str(index/obj.fy * 1e3) ' mm']);
                    xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
                    ylabel('\fontname {times new roman} z (mm)', 'fontsize', 16);
                    mask  = (obj.row==index) & obj.col>=win(1) & obj.col<=win(end);
                    x_row = obj.col(mask) / obj.fx * 1e3;
                    x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000/2;  % / obj.fy * 1e3 ;
                    h1 = scatter(ca1, x_row, x_dep, ...
                        3,  [122 122 121]/255, 'filled', ...
                        'DisplayName','Interply track');
                    hold on;
                    h2 = scatter(x, front_p(index, win)/ obj.fs * 1e3 * 3000/2, ...
                        3, 'red', 'filled', ...
                        'DisplayName','Front surface');
                    hold on;
                    h3 = scatter(x, bakc_p(index, win)/ obj.fs * 1e3 * 3000/2, ...
                        3, 'magenta', 'filled', ...
                        'DisplayName','Rear surface');
                    hold on;
                    % figure setup
                    colormap(hsv);
                    caxis(ca1, [-90 90]);
                    h = colorbar;
                    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (\circ)');
                    set(ca1, 'Fontname', 'times new Roman');
                    set(ca1, 'Fontsize', 16);
                    set(ca1, 'Linewidth', 1.5);
                    % top view
                    figure(cf);
                    ca2     = subplot(2, 1, 2);
                    % imagesc(x, z, B_scan');
                    z_index     = round(mean(front_p, 'all', 'omitnan') ...
                        + mean(bakc_p-front_p, 'all', 'omitnan') / size(dataset,1) * index);
                    B_scan = squeeze(dataset(1:end, 1:end, z_index));
                    pcolor(x, y, B_scan);
                    shading flat;
                    hold on;
                    title(['z =' num2str(z_index/obj.fs * 1e3 * 3000/2) ' mm']);
                    xlabel('\fontname {times new roman} x (mm)', 'fontsize', 16);
                    ylabel('\fontname {times new roman} y (mm)', 'fontsize', 16);
                    hold on;
                    % figure setup
                    colormap(hsv);
                    caxis(ca2, [-90 90]);
                    h = colorbar;
                    set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (\circ)');
                    set(ca2, 'Fontname', 'times new Roman');
                    set(ca2, 'Fontsize', 16);
                    set(ca2, 'Linewidth', 1.5);
                    % frame vedio operation
                    frame = getframe(cf);
                    writeVideo(v, frame);
                    clf(cf);
                end
            elseif (B_type == 'y')
            end
        end
        
        % ************************ Global orientation extract ********************
        function obj = extract_gloabl_orientation_3D(obj, PropertyName, wavelength, orientation, ratios, max_layer, SFB, SAR, K)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ratios: the depth in each ply
            % max_layer: the total plies of the structure
            % SFB, SAR: SpatialFrequencyBandwidth and SpatialAspectRatio
            % K: the smoothing applied to the Gabor magnitude responses.
            % ***
            gaborArray             = gabor(wavelength, orientation, ...
                'SpatialFrequencyBandwidth', SFB, 'SpatialAspectRatio', SAR);
            wl_max                 = max(wavelength);
            Inplane_direction      = NaN(size(obj.(PropertyName)));
            Inp_direc_domain       = NaN(size(obj.(PropertyName)));
            Inp_direc_domain_layer = NaN(size(obj.(PropertyName)));
            ratio_spacing          = (ratios(end) - ratios(1)) / (length(ratios) - 1);
            for i = 1:max_layer
                ID_sum_layer = zeros(length(wavelength), length(orientation));
                for j = 1:length(ratios)
                    ratio                                          = ratios(j);
                    [Cscan_inam, profile_layer, profile_layer_pre] = obj.define_Cscan_depthprofile_knn(i, ratio, PropertyName, max_layer);
                    %                     % fillna and smooth the profile
                    %                     profile_layer                                  = fx_inpaint_nans(profile_layer, 5);
                    %                     profile_layer_pre                              = fx_inpaint_nans(profile_layer_pre, 5);
                    profile_layer                                  = medfilt2(profile_layer, [7, 7]);
                    profile_layer_pre                              = medfilt2(profile_layer_pre, [7, 7]);

                end
                % end loop in one Cscan
                disp([num2str(i), '/', num2str(max_layer), ' layer']);
            end
            obj.Inplane_direction_3D               = Inplane_direction;
            obj.Inplane_direction_3D_overall       = Inp_direc_domain;
            obj.Inplane_direction_3D_overall_layer = Inp_direc_domain_layer;
        end
         
        % *********************** kNN clustering for interply track *********************
        function obj = show_interply_track_3D(obj)
            % convert the mask_plytrack to ordinates of 3d scattering points
            % display the 3d scattering with the layer marks
            % initial the props for clustering
            %             temp = obj.mask_plytrack;
            %             obj.row = [];
            %             obj.col = [];
            %             obj.dep = [];
            %             obj.layer_counts = [];
            %             for i = 1:size(temp, 1)
            %                 for j = 1:size(temp, 2)
            %                     % find the sequence of the '1' in the mask along the depth dimension.
            %                     idx = find(temp(i, j, :));
            %                     obj.layer_counts = cat(2, obj.layer_counts, (1:length(idx)));
            %                     obj.row = cat(2, obj.row, i*ones(1, length(idx)));
            %                     obj.col = cat(2, obj.col, j*ones(1, length(idx)));
            %                     obj.dep = cat(2, obj.dep, idx');
            %                 end
            %                 disp([num2str(i) '/' num2str(size(temp, 1))]);
            %             end
            cf = figure();
            scatter3(obj.row, obj.col, obj.dep, 10, obj.layer_counts.^4, 'o', 'filled');
            set(gca, 'ZDir', 'reverse');
            obj.X_scatter = [obj.row; obj.col; obj.dep]';
            zlim([1000 1250]);
        end
        
        function obj = knn_search(obj, K, P, nbins)
            % search k nearest neighbors and display the hist of distance
            % K: number of the neighbors
            % parameters of the Minkowski distance
            % nbins: number of bins for the hist
            [obj.cIdx, obj.cD] = knnsearch(obj.X_scatter, obj.X_scatter, 'K', K, 'NSMethod', 'kdtree', 'Distance', 'Minkowski', 'P', P);
            % hist the distances
            cf = figure();
            h = histogram(mean(obj.cD, 2), nbins);
            set(gca,'YScale','log');
            xlabel('\fontname {times new roman} Average distance to k nearest neighbors (arb.) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Count', 'fontsize', 16);
        end
        
        function obj = knn_clustering(obj, threshold_dis)
            % clustering the interply track by knn
            % threshold_dis: the distance threshold to determine the noisy points
            % clustering & remove the noisy points
            obj.layer_counts_new = obj.layer_counts;
            for i = 1:length(obj.layer_counts)
                dis_ave = mean(obj.cD(i, :));
                if dis_ave > threshold_dis
                    obj.layer_counts_new(i) = NaN;
                else
                    obj.layer_counts_new(i) = mode(obj.layer_counts(obj.cIdx(i, :)));
                end
                if mod(i,10000)==0
                    disp([num2str(i) '/' num2str(length(obj.layer_counts))]);
                end
            end
            cf = figure();
            scatter3(obj.row, obj.col, obj.dep, 10,  obj.layer_counts_new.^4, 'o', 'filled');
            set(gca, 'ZDir', 'reverse');
            zlim([1000 1250]);
        end
        
        function obj = clear_knn_properties(obj)
            % clear the props for KNN to release the memory
            obj.layer_counts = [];
            obj.X_scatter = [];
            obj.cIdx = [];
            obj.cD = [];
        end
        
        function show_one_interply_knn(obj, layer, v)
            % display one interply track in C-scan mode after knn
            % layer: 'layer'th layer
            % v: the average velocity of the laminate
            oneply_row = obj.row(obj.layer_counts_new==layer);
            oneply_col = obj.col(obj.layer_counts_new==layer);  % / obj.fy * 1e3 ;
            oneply_dep = obj.dep(obj.layer_counts_new==layer);  % / obj.fy * 1e3 ;
            oneply_dep_2D = NaN(max(oneply_row), max(oneply_col));
            for i = 1: length(oneply_row)
                oneply_dep_2D(oneply_row(i), oneply_col(i)) = (oneply_dep(i) - ...
                    obj.front_I(oneply_row(i), oneply_col(i))) / obj.fs * 1e3 * v / 2;
            end
            cf = figure();
            ax = subplot(1, 1, 1);
            x = (1:size(oneply_dep_2D, 1)) / obj.fx * 1e3;
            y = (1:size(oneply_dep_2D, 2)) / obj.fy * 1e3;
            imagesc(ax, x, y, oneply_dep_2D);
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Depth (mm)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function obj = show_interply_track_3D_knn(obj, xslice, yslice, zslice)
            % demonstrate the 3d ply track after knn
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % ***
            inph_ex = abs(obj.img_hil);
            y = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e3 * 3000 / 2;
            [X, Y, Z] = meshgrid(x, y, z);
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure('Name', ['3dtrack_interply_' num2str(xslice(1)*obj.fx/1e3)]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            h = slice(ax, X, Y, Z, inph_ex, xslice, yslice, zslice);
            set(get(get(h(1),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            set(get(get(h(2),'Annotation'),'LegendInformation'),'IconDisplayStyle','off');
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap gray;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Z (mm)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            xlim([0 x(end)]);
            ylim([0 y(end)]);
%             zlim([0.8 7]);
            zticks([1 3 5]);
            zticklabels({'0','2','4'});
            % select the xslice
            x_idx = round(xslice * obj.fx / 1e3);
            y_idx = round(yslice * obj.fy / 1e3);
            mask = (obj.col==x_idx) | (obj.row==y_idx);
            x_row = obj.row(mask) / obj.fx * 1e3;
            x_col = obj.col(mask) / obj.fy * 1e3;  % / obj.fy * 1e3 ;
            x_dep = obj.dep(mask) / obj.fs * 1e3 * 3000 / 2;  % / obj.fy * 1e3 ;
            scatter3(ax, x_col, x_row, x_dep, ...
                2, 'cyan', 'filled', ...
                'DisplayName','Interply track');
            hold on;
            scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                y, obj.front_I(:, x_idx)/ obj.fs * 1e3 * 3000 / 2, ...
                3, 'red', 'filled', ...
                'DisplayName','Front surface');
            hold on;
            scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                y, obj.rear_I(:, x_idx)/ obj.fs * 1e3 * 3000 / 2, ...
                3, 'magenta', 'filled', ...
                'DisplayName','Rear surface');
            hold on;
            % select the yslice
            scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                obj.front_I(y_idx, :)/ obj.fs * 1e3 * 3000 / 2, ...
                3, 'red', 'HandleVisibility','off');
            hold on;
            scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                obj.rear_I(y_idx, :)/ obj.fs * 1e3 * 3000 / 2, ...
                3, 'magenta', 'HandleVisibility','off');
            hold off;
            legend;
            view([15 65 40]);
        end
        
        % finish tne methods
    end
    
end


