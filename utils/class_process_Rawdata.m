classdef class_process_Rawdata
    % Addme
    % process the cleaned RF raw volumeric data
    % a python script has been recommended to preprocess the data
    % then a cleaned data are produced including:
    % img_pre: preprocessed img
    % img_hil: analytical signal of preprocessed img
    % front: preprocessed front surface index
    % rear: preprocessed rear surface index
    
    properties
        filename
        img
        img_hil
        Props
        setup
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
        % the object of reference signal: 'class_reference_signal'
        reference
        refer_Ascan_ori
        refer_Ascan
        refer_Ascan_aligned
        % deconvolution
        img_WienerDeconv_B
        img_WienerDeconv_AR_B
        img_SparseDeconv_B
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
        RT_ID
        RT_ID_gyy
        rear_mask
        Inplane_direction_3D
        Inplane_direction_3D_overall
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
    end
    
    methods
        % ***************** definition, load data **************
        function obj = class_process_Rawdata(filename)
            % filename: path and name of the mat. data including: img_pre, img_hil, front, rear
            obj.filename = filename;
        end
        
        function obj = read_data(obj)
            my_tdms_struct             = TDMS_getStruct(obj.filename);
            [obj.Props, obj.setup]   = fx_read_setup(my_tdms_struct);
            
            disp ('start reading');
            % read raw data from tdms_struct
            index_cell = fieldnames(my_tdms_struct);
            for i = 2:numel(fieldnames(my_tdms_struct))
                fieldname = index_cell{i};
                index_struct = my_tdms_struct.(fieldname);
                scan_point_cell = fieldnames(index_struct);
                for j = 3:numel(fieldnames(index_struct))
                    fieldname = scan_point_cell{j};
                    scan_point_struct = index_struct.(fieldname);
                    if (i==2)&&(j==3)
                        oimg_uint8 = zeros(numel(fieldnames(my_tdms_struct))-1,...
                            numel(fieldnames(index_struct))-1, ...
                            length(scan_point_struct.data));  
                        oimg_uint8(i-1, j-2, :) = scan_point_struct.data;
                    else
                        oimg_uint8(i-1, j-2, :) = scan_point_struct.data;
                    end
                end
            end
            obj.img           = oimg_uint8;
            img_hil_temp = oimg_uint8;
            disp ('start hilbert transform');
            for i = 1:size(oimg_uint8, 1)
                for j = 1:size(oimg_uint8, 2)
                    img_hil_temp(i, j, :) = hilbert(oimg_uint8(i, j, :));
                end
            end
             obj.img_hil = img_hil_temp;
            obj.fs           = str2double(obj.setup.SampleRateHz); % read the samping frequency
            obj.fx           = 1e3 / str2double(obj.setup.Breakpointscanmm); % read the x step and turn it into fx: mm -> 1 / m
            obj.fy           = 1e3 / str2double(obj.setup.Breakpointindexmm); % read the y step and turn it into fy: mm -> 1 / m
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
        
        % ************************************
        function show_hilbert_Ascan(obj, x, y)
            % display the decomposed analytical signal of A scan
            % x, y: index of the A scan
            % reuse the fx_showAS
            ascan = squeeze(obj.img(x, y, :));
            ascan_hil = hilbert(ascan);
            fx_showAS((1:length(ascan)) / obj.fs, ascan_hil, obj.fs);
        end
        
        function show_Ascan_inam_peaks(obj, x, y, MinPD, MinPH)
            % demonstrate the inam of a ascan
            % find the peaks in the inam
            % x: x index
            % y: y index
            % MinPD: MinPeakDistance for findpeaks
            % MinPH: MinPeakHeight for findpeaks
            % ***
            % reset front surfaces
            obj.front_I = obj.front;
            % calculate the inam
            % The sequence is changed unexpectedly. Transverse the matrix
            inam = abs(obj.img_hil);
            inph = angle(obj.img_hil);
            %
            t = (1:size(inam, 3));
            A_inam = squeeze(inam(x, y, :));
            A_inph = squeeze(inph(x, y, :));
            figure;
            ca = subplot(2, 1, 1);
            plot(t, A_inam, 'linewidth', 2);
            %             xlabel('\fontname {times new roman} Time(\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Inst. amp.', 'fontsize', 16);
            hold on;
            % find peaks
            [pks, locs, ~, ~] = findpeaks(A_inam, t, 'MinPeakDistance', MinPD, 'MinPeakHeight', MinPH);
            % select 2 maximum peaks
            [~, I] = sort(pks);
            locs_max = locs(I);
            plot(t(locs_max), A_inam(locs_max), 'rv');
            hold on;
            ca = subplot(2, 1, 2);
            plot(t, A_inph, 'linewidth', 2);
            %             xlabel('\fontname {times new roman} Time(\mus) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Inst. amp.', 'fontsize', 16);
            hold on;
            plot(t(locs_max), A_inph(locs_max), 'rv');
            hold on;
        end
        
        function show_Cscan(obj, z, PropertyName)
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
            imagesc(ax, X,Y, C_scan / obj.fx * 1e6)
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            %
            figure('Name', ['C_scan_inph' , '_', num2str(z)]);
            C_scan = angle(obj.(PropertyName)(:, :, z));
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan, 2)) / obj.fy * 1e3;
            imagesc(ax, X,Y, C_scan / obj.fx * 1e6)
            hold on;
            h = colorbar;
            %             caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Time(\mus)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
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
        
        function demo_Ascan(obj, x, y, fig_subname)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % fig_subname: the subname for the fig
            Ascan = squeeze(obj.img(x, y, :));
            % remove the direct component 
            Ascan = Ascan - mean(Ascan);
            % plot the inam and the record signal
            cf = figure('Name', strcat('Ascan_x', num2str(x), '_y', num2str(y), fig_subname) );
            set(cf, 'Position', [0, 0, 1000, 600], 'color', 'white');
            ca = subplot(1, 1, 1);
            plot(ca, (1:length(Ascan)) / obj.fs * 1e6, Ascan, 'linewidth', 2);
            hold on;
            plot(ca, (1:length(Ascan)) / obj.fs * 1e6, abs(hilbert(Ascan)), 'linewidth', 2);
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
        
        % ***************** surface determination ****************
        function obj = find_front_amp(obj,  MinPD, MinPH, PropertyName)
            % determining the front surface index by amp.
            % obj.rear_i is also intially determined here
            % The sequence is changed unexpectedly. Transverse the matrix
            % MinPD: MinPeakDistance for findpeaks
            % MinPH: MinPeakHeight for findpeaks
            % PropertyName: 'img_hil' or 'img_hil_filter' ... used for the inam.
            % ****
            % calculate the inam
            % The sequence is changed unexpectedly. Transverse the matrix
            inam = abs(obj.(PropertyName));
            inph = angle(obj.(PropertyName));
            % the time domain
            t = (1:size(inam, 3));
            front_I_temp = NaN(size(inam, 1), size(inam, 2));
            rear_I_temp = NaN(size(inam, 1), size(inam, 2));
            %
            for i = 1:size(inam, 1)
                for j = 1:size(inam, 2)
                    A_inam = squeeze(inam(i, j, :));
                    A_inph = squeeze(inph(i, j, :));
                    [pks, locs, ~, ~] = findpeaks(A_inam, t, 'MinPeakDistance', MinPD, 'MinPeakHeight', MinPH);
                    %                     % sort the peaks
                    %                     [~, I] = sort(pks);
                    % the front-wall echo is largest
                    if length(pks)>1
                        % the first echo should be front
                        front_I_temp(i, j) = min(locs);
                        % derive the front_I by phase
                        A_inph_locs = abs(A_inph(locs) - A_inph(min(locs)));
                        [~, I_phase] = sort(A_inph_locs);
                        % find the fisrt element larger than the front echo
                        if find(pks>pks(1), 1)
                            rear_I_temp(i, j) = locs(find(pks>pks(1), 1));
                        elseif locs(I_phase(1)) < locs(I_phase(end))
                            % front_I_temp(i, j) = locs(I_phase(1));
                            rear_I_temp(i, j) = locs(I_phase(end));
                        end
                    elseif length(pks)==1
                        locs_max = locs(1);
                        front_I_temp(i, j) = locs_max;
                    end
                end
                disp([num2str(i) '/' num2str(size(inam, 1))]);
            end
            obj.front_I = front_I_temp;
            obj.rear_I = rear_I_temp;
        end
        
        function obj = smooth_rear_I(obj, ~)
            % fill and smooth the rear_I
            obj.front_I_pre = obj.front_I;
            obj.rear_I_pre = obj.rear_I;
            obj.front_I = fx_inpaint_nans(obj.front_I, 5);
            obj.front_I = round(medfilt2(obj.front_I, [5, 5], 'symmetric'));
            obj.rear_I = fx_inpaint_nans(obj.rear_I, 5);
            obj.rear_I = round(medfilt2(obj.rear_I, [5, 5], 'symmetric'));
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
            rear_I_nan = obj.rear_I;
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
        
        function show_surfaces(obj)
            % display the surfaces: phase, amp., index
            % calculate the inam and inph
            % The sequence is changed unexpectedly. Transverse the matrix
            inph = atan2(imag(obj.img_hil), real(obj.img_hil));
            %             inph = permute(inph, [2 3 1]);
            inam = abs(obj.img_hil);
            %             inam = permute(inam, [2 3 1]);
            % phase calculation
            phase_front = zeros(size(obj.front_I));
            amp_front = zeros(size(obj.front_I));
            phase_rear = zeros(size(obj.front_I));
            amp_rear = zeros(size(obj.front_I));
            %             phase_front_sum = 0;
            front_I_mean = round(mean(obj.front_I, 'all', 'omitnan'));
            rear_I_mean = round(mean(obj.rear_I, 'all', 'omitnan'));
            %             rear_I_mean = round(mean(mean(obj.rear_I)));
            for i = 1:size(inam, 1)
                for j = 1:size(inam, 2)
                    % if front index is out of boundary, set mean
                    if isnan(obj.front_I(i, j)) || round(obj.front_I(i, j)) <= 0 || round(obj.front_I(i, j)) > size(inph, 3)
                        phase_front(i, j) = inph(i, j, front_I_mean);
                        amp_front(i, j) = inam(i, j, front_I_mean);
                    else
                        phase_front(i, j) = inph(i, j, round(obj.front_I(i, j)));
                        amp_front(i, j) = inam(i, j, round(obj.front_I(i, j)));
                    end
                    % if rear index is out of boundary, set mean
                    if isnan(obj.rear_I(i, j)) || round(obj.rear_I(i, j)) > size(inph, 3) || round(obj.rear_I(i, j)) <= 0
                        phase_rear(i, j) = inph(i, j, rear_I_mean);
                        amp_rear(i, j) = inam(i, j, rear_I_mean);
                    else
                        phase_rear(i, j) = inph(i, j, round(obj.rear_I(i, j)));
                        amp_rear(i, j) = inam(i, j, round(obj.rear_I(i, j)));
                    end
                end
            end
            %             phase_rear_mean = phase_front_mean + pi; % rad
            % visual the phase_front amp_front
            figure, ca = subplot(3, 1, 1);
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
        
        function show_inam_3d(obj, PropertyName, xslice, yslice, zslice)
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
            inph             = abs(obj.(PropertyName));
            inph_ex       = inph;
            y                  = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x                  = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z                  = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z]     = meshgrid(x, y, z);
            inph_visual = inph_ex;
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            ax                = subplot(1, 1, 1);
            h                  = slice(ax, X, Y, Z, inph_visual, xslice, yslice, zslice);
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h                   = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Inst. phase(rad.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            set(gca,'zDir','reverse'); % get time domain upside down
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
        end
        
        function obj = cut_edges(obj, window)
            % cut the edges of the dataset, and display the processed
            % figure
            % window: [x1, x2, y1, y2, t1, t2]
            x1 = window(1);
            x2 = window(2);
            y1 = window(3);
            y2 = window(4);
            t1 = window(5);
            t2 = window(6);
            obj.img = obj.img(x1: min(end, x2), y1: min(end, y2), t1: min(end, t2));
            obj.img_hil = obj.img_hil(x1: min(end, x2), y1: min(end, y2), t1: min(end, t2));
            obj.front = obj.front(x1:min(end, x2), y1:min(end, y2)) - t1;
            obj.rear = obj.rear(x1:min(end, x2), y1:min(end, y2)) - t1;
            obj.front_I = obj.front_I(x1:min(end, x2), y1:min(end, y2)) - t1;
            obj.rear_I = obj.rear_I(x1:min(end, x2), y1:min(end, y2)) - t1;
        end
        
        function obj = perfect_rear_surface(obj)
            % make the rear_surface to the end of the 3D volumetric dataset
            obj.rear_I = size(obj.img, 3) * ones(size(obj.img, 1), size(obj.img, 2));
            obj.rear_mask = [];
        end
        
        % ******************** filter ***********
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
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % track the - pi / 2 by interplation on the unwraped data
            A_inph = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            A_inph_unw = unwrap(A_inph);
            phase_seq = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
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
            % threshold: the threshold to select the position of inst. phases
            Ascan          = squeeze(obj.img(x, y, :));
            front_I_interp = obj.front_I(x, y);
            rear_I_interp  = obj.rear_I(x, y);
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
            filter         =  fx_1dLogGabor_filter(fshift, f0, sigma);
            Fa_filter      = Fa .* filter;
            % ifft to obtain the AS
            Fa_ishift      = ifftshift(Fa_filter);  % back to the original transform output
            Ascan_as       = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as       = conj(flip(Ascan_as));
            % plot the inam and the record signal
            cf             = figure('Name', strcat('plytrack_filter', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca             = subplot(2, 1, 1);
            t_space        = (1:length(Ascan_as)) / obj.fs;
            plot(ca, t_space, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, t_space, abs(Ascan_as), 'linewidth', 2)
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Recorded signal', 'Inst. Amp.'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % track the - pi / 2 by interplation on the unwraped data
            A_inph                                  = angle(Ascan_as);
            %             A_inph_win                            = A_inph(round(front_I_interp): round(rear_I_interp)); % inph in the surface window
            I_phase_track                           = abs((A_inph + pi / 2)) < threshold;
            I_phase_track(1:round(front_I_interp))  = 0;
            I_phase_track(round(rear_I_interp):end) = 0;
            % plot the inph, and the ply track of - pi / 2
            ca = subplot(2, 1, 2);
            plot(ca, t_space, A_inph, 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, [front_I_interp rear_I_interp] / obj.fs , angle(Ascan_as(round([front_I_interp rear_I_interp]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            hold on;
            plot(ca, t_space(I_phase_track), A_inph(I_phase_track) , 'gd', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Inst. phase', 'Derived surface', 'Derived inter-ply'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
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
            Ascan_as = ifft(Fa_ishift);
            % the direction of 'as' is reversed somehow, needs to be
            % flipped and conjugated!
            Ascan_as = conj(flip(Ascan_as));
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
            A_inph = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
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
            Ascan = squeeze(obj.img(x, y, :));
            front_I_interp = obj.front_I(x, y);
            rear_I_interp = obj.rear_I(x, y);
            % add noise
            Ascan = awgn(Ascan, snr_value, 'measured');
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
            % calcualte the snr
            %             SNR = snr(real(Ascan_as), obj.fs);
            pn = bandpower(real(Ascan_as), obj.fs, [10e6, obj.fs / 2.1]);
            ps = bandpower(real(Ascan_as), obj.fs, [1e6, 10e6]);
            SNR = 10 * log10(ps / pn);
            disp("SNR:");
            disp(SNR);
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
            legend('Location', 'best');
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % track the - pi / 2 by interplation on the unwraped data
            A_inph = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            A_inph_unw = unwrap(A_inph);
            phase_seq = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
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
            legend('Location', 'northwest');
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        % ***************** structural tensor to extract angle ****************
        function obj = structural_tensor(obj, sigma1, sigma2, PropertyName)
            % apply structural tensor on the cleaned inph dataset
            % explanation about the angles:
            % anglex, y, z represent the angles with the axis directional vector of the planar like structure.
            % extract the angles of the the plane-structrual
            % sigma1: smoothing scale
            % sigma2: integration scale
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % ***
            % time record
            tic;
            disp('Start structural tensor method');
            % The sequence is changed unexpectedly. Transverse the matrix
            inph = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            % 3D Gaussian filter
            inph_cos = cos(inph);
            inph_sin = sin(inph);
            inph_cos_filter = imgaussfilt3(inph_cos, sigma1, 'padding', 'replicate');
            inph_sin_filter = imgaussfilt3(inph_sin, sigma1, 'padding', 'replicate');
            % calculation of gradients
            [Gx_cos, Gy_cos, Gz_cos] = fx_imgradientxyz(inph_cos_filter, 'intermediate');
            [Gx_sin, Gy_sin, Gz_sin] = fx_imgradientxyz(inph_sin_filter, 'intermediate');
            Gx = inph_cos.*Gx_sin - inph_sin.*Gx_cos;
            Gy = inph_cos.*Gy_sin - inph_sin.*Gy_cos;
            Gz = inph_cos.*Gz_sin - inph_sin.*Gz_cos;
            % release the memory
            clear inph_cos;
            clear inph_sin;
            clear inph;
            clear inam;
            clear infq;
            clear inph_cos_filter;
            clear inph_sin_filter;
            % for debug
            % [~, ~] = fx_showBS(2, inph, num_bs, step_bs, data.front, data.rear, B_type, 'jet');
            %             [~, ~] = fx_showBS(3, Gx, num_bs, step_bs, data.front, data.rear, B_type, 'jet');
            %             [~, ~] = fx_showBS(4, Gy, num_bs, step_bs, data.front, data.rear, B_type, 'jet');
            %             [fig, ax] = fx_showBS(5, Gz, num_bs, step_bs, data.front, data.rear, B_type, 'jet');
            % form the structure-tensor
            ST = fx_formStructureTensor(Gx, Gy, Gz);
            ST_in = zeros(size(ST));
            for i=1:3
                for j=1:3
                    ST_in(:, :, :, i, j) = imgaussfilt3(ST(:, :, :, i, j), sigma2);
                end
            end
            % time record
            timeElapsed = toc;
            disp(['form_tensor: ', num2str(timeElapsed)]);
            [~, cp, anglex, angley, anglez] = fx_decomSS_surface(ST_in, obj.front_I, obj.rear_I);
            % release the memory
            clear Gx;
            clear Gy;
            clear img;
            clear img_deconv;
            clear data;
            clear inam;
            % time record
            timeElapsed = toc;
            disp(['extraxt angles: ', num2str(timeElapsed)]);
            obj.c_p = cell2mat(cp);
            obj.angle_x = cell2mat(anglex);
            obj.angle_y = cell2mat(angley);
            obj.angle_z = cell2mat(anglez);
            % time record
            timeElapsed = toc;
            disp(['convert to arrays: ', num2str(timeElapsed)]);
        end
        
        function show_angles_ST(obj, medf_kernel, xslice, yslice, zslice)
            % display the extracted angles of ST
            % filter the noise by medfilter
            % medf_kernel: the kernel of median filter for smoothing
            % xslice, yslice, zslice: the slice in the 3d visualization;
            % these should be real unit here: mm, mm, us
            % for example:   xslice =  2; yslice = [2 5]; zslice = [];
            cp = medfilt3(obj.c_p, medf_kernel, 'replicate');
            anglex = medfilt3(real(obj.angle_x), medf_kernel, 'replicate');
            angley = medfilt3(real(obj.angle_y), medf_kernel, 'replicate');
            % anglex
            y = (0: size(anglex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(anglex, 2) - 1) / obj.fy * 1e3;
            z = (size(anglex, 3) - 1: -1: 0) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x,y,z);
            cf = figure; h = slice(X, Y, Z, anglex, xslice, yslice, zslice);
            ax = subplot(1, 1, 1);
            set(h,'EdgeColor','none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle(degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            title('X-Z angle', 'fontsize', 16, 'Fontname', 'times new Roman');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            % angley
            cf = figure; h = slice(X, Y, Z, angley, xslice, yslice, zslice);
            ax = subplot(1, 1, 1);
            set(h,'EdgeColor','none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle(degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            title('Y-Z angle', 'fontsize', 16, 'Fontname', 'times new Roman');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            % c_p
            cf = figure; h = slice(X, Y, Z, cp, xslice, yslice, zslice);
            ax = subplot(1, 1, 1);
            set(h,'EdgeColor','none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle(degree)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            zlabel('\fontname {times new roman} Time (\mus)', 'fontsize', 16);
            title('C_p', 'fontsize', 16, 'Fontname', 'times new Roman');
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            % caxis([-0.5, 0.5]);
        end
        
        function obj = Filter_logGabor(obj, f0, sigma, PropertyName)
            % Filter the 3d dataset by logGabor filter
            % f0: the array of the center frequencies of the filter
            % sigma: affects the bandwidth of the filter.
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            [lx, ly, n] = size(obj.(PropertyName));
            % save in the obj.img_hil_filter
            obj.img_hil_filter = zeros(size(obj.(PropertyName)));
            % define a logGabor filter and calcualte the analytical signal
            fshift = (-n/2: n/2-1) * (obj.fs / n); % zero-centered frequency range
            filter =  fx_1dLogGabor_filter(fshift, f0, sigma);
            for x = 1: lx
                for y = 1:ly
                    Ascan = squeeze(obj.(PropertyName)(x, y, :));
                    % calculate the AS in frequency domain
                    Fo = fft(Ascan);
                    % transverse Fo for dot multiply
                    Fo = Fo';
                    fx_shift = fftshift(Fo);
                    Fa = (1 + sign(fshift)) .* fx_shift;
                    % add the logGabor filter here
                    Fa_filter = Fa .* filter;
                    % ifft to obtain the AS
                    Fa_ishift = ifftshift(Fa_filter);  % back to the original transform output
                    Ascan_as = ifft(Fa_ishift);
                    % the direction of 'as' is reversed somehow, needs to be
                    % flipped and conjugated!
                    Ascan_as = conj(flip(Ascan_as));
                    obj.img_hil_filter(x, y, :) = Ascan_as;
                end
            end
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
            inph = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            mask_plytrack_temp = zeros(size(inph));
            row_temp = [];
            col_temp = [];
            dep_temp = [];
            layer_counts_temp = [];
            for idx = 1:size(inph, 1)
                for idy = 1:size(inph, 2)
                    if obj.rear_I(idx, idy) <=0 || obj.rear_I(idx, idy) > size(inph, 3) ...
                            || obj.front_I(idx, idy) <=0 || obj.front_I(idx, idy) > obj.rear_I(idx, idy)
                        continue;
                    end
                    A_inph = squeeze(inph(idx, idy, round(obj.front_I(idx, idy)): ...
                        round(obj.rear_I(idx, idy))));
                    % track the phase by unwraped inst. phase
                    A_inph_unw = unwrap(A_inph);
                    phase_seq = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
                    I_phase_track = interp1(A_inph_unw, 1:length(A_inph_unw), phase_seq)...
                        + round(obj.front_I(idx, idy));
                    mask_plytrack_temp(idx, idy, round(I_phase_track)) = 1;
                    % derive the rear surface by phase
                    A_inph = squeeze(inph(idx, idy, round(obj.front_I(idx, idy)):end));
                    A_inph_unw = unwrap(A_inph);
                    phase_rear_surface = phase_seq(end) + 3 * pi / 2;  % -pi for rear surface, which means 3 pi / 2 from the last interply
                    I_rear_surface = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                        phase_rear_surface) + round(obj.front_I(idx, idy));
                    obj.rear_I(idx, idy) = round(I_rear_surface);
                    % save the 3D scatter points
                    layer_counts_temp = cat(2, layer_counts_temp, (1:length(I_phase_track)));
                    row_temp = cat(2, row_temp, idx * ones(1, length(I_phase_track)));
                    col_temp = cat(2, col_temp, idy * ones(1, length(I_phase_track)));
                    dep_temp = cat(2, dep_temp, round(I_phase_track));
                end
                disp([num2str(idx), '/' , num2str(size(inph, 1))]);
            end
            obj.mask_plytrack = mask_plytrack_temp;
            obj.row = row_temp;
            obj.col = col_temp;
            obj.dep = dep_temp;
            obj.layer_counts = layer_counts_temp;
        end
        
        function obj = track_interply_inph(obj, threshold, PropertyName)
            % track the interply by instantaneous phases.
            % also derive the rear_I by phase
            % return a masked dataset
            % threshold: the threshold to select the position of inst. phases
            % phase close to -pi/2
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % ************
            % The sequence is changed unexpectedly. Transverse the matrix
            inph               = atan2(imag(obj.(PropertyName)), real(obj.(PropertyName)));
            mask_plytrack_temp = zeros(size(inph));
            %             row_temp           = [];
            %             col_temp           = [];
            %             dep_temp           = [];
            %             layer_counts_temp  = [];
            for idx = 1:size(inph, 1)
                for idy = 1:size(inph, 2)
                    if obj.rear_I(idx, idy) <=0 || obj.rear_I(idx, idy) > size(inph, 3) ...
                            || obj.front_I(idx, idy) <=0 || obj.front_I(idx, idy) > obj.rear_I(idx, idy)
                        continue;
                    end
                    A_inph                            = squeeze(inph(idx, idy, round(obj.front_I(idx, idy)): ...
                        round(obj.rear_I(idx, idy))));
                    % track the phase by wrapped inst. phase
                    I_phase_track                     = abs(A_inph + pi/2) < threshold;
                    mask_plytrack_temp(idx, idy, round(obj.front_I(idx, idy):...
                        round(obj.rear_I(idx, idy)))) = I_phase_track;
                    % derive the rear surface by phase
                    I_rear_surface                    = obj.rear_I(idx, idy);  % -pi for rear surface
                    %                     I_rear_surface       = interp1(A_inph_unw, 1:length(A_inph_unw), ...
                    %                         phase_rear_surface) + round(obj.front_I(idx, idy));
                    obj.rear_I(idx, idy)              = round(I_rear_surface);
                    %                     % save the 3D scatter points
                    %                     layer_counts_temp = cat(2, layer_counts_temp, (1:length(I_phase_track)));
                    %                     row_temp = cat(2, row_temp, idx * ones(1, length(I_phase_track)));
                    %                     col_temp = cat(2, col_temp, idy * ones(1, length(I_phase_track)));
                    %                     dep_temp = cat(2, dep_temp, round(I_phase_track));
                end
                disp([num2str(idx), '/' , num2str(size(inph, 1))]);
            end
            obj.mask_plytrack = mask_plytrack_temp;
%             obj.row = row_temp;
%             obj.col = col_temp;
%             obj.dep = dep_temp;
%             obj.layer_counts = layer_counts_temp;
        end
        
        function demo_Ascan_plytrack(obj, x, y, interp_factor)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % interp_factor: the interpolation multiple
            Ascan = squeeze(obj.img(x, y, :));
            % interpolation
            xp = 1: 1: length(Ascan);
            xq = 1: 1 / interp_factor: length(Ascan);
            Ascan_interp = interp1(xp, Ascan, xq);
            front_I_interp = obj.front_I(x, y) * interp_factor;
            rear_I_interp = obj.rear_I(x, y) * interp_factor;
            % plot the inam and the record signal
            Ascan_as = hilbert(Ascan_interp);
            cf = figure('Name', strcat('plytrack', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(2, 1, 1);
            plot(ca, (1:length(Ascan_interp)) / obj.fs, Ascan_interp, 'linewidth', 2);
            hold on;
            plot(ca, (1:length(Ascan_as)) / obj.fs, abs(Ascan_as), 'linewidth', 2)
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Recorded', 'Inst. Amp.'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % track the - pi / 2 by interplation on the unwraped data
            A_inph = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            A_inph_unw = unwrap(A_inph);
            phase_seq = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
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
        
        function demo_Ascan_plytrack_inph(obj, x, y, interp_factor, threshold)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % interp_factor: the interpolation multiple
            % 
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
            t_space        = (1:length(Ascan_as)) / obj.fs;
            plot(ca, t_space, real(Ascan_as), 'linewidth', 2);
            hold on;
            plot(ca, t_space, abs(Ascan_as), 'linewidth', 2)
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Recorded signal', 'Inst. Amp.'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % track the - pi / 2 by interplation on the unwraped data
            A_inph                                  = angle(Ascan_as);
%             A_inph_win                            = A_inph(round(front_I_interp): round(rear_I_interp)); % inph in the surface window
            I_phase_track                           = abs((A_inph + pi / 2)) < threshold;
            I_phase_track(1:round(front_I_interp))  = 0;
            I_phase_track(round(rear_I_interp):end) = 0;
            % plot the inph, and the ply track of - pi / 2
            ca = subplot(2, 1, 2);
            plot(ca, t_space, A_inph, 'linewidth', 2); % angle(Ascan_as) is the full inph
            hold on;
            plot(ca, [front_I_interp rear_I_interp] / obj.fs , angle(Ascan_as(round([front_I_interp rear_I_interp]))) ...
                , 'ro', 'linewidth', 2); % the surfaces
            hold on;
            plot(ca, t_space(I_phase_track), A_inph(I_phase_track) , 'gd', 'linewidth', 2);
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Phase (rad.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Inst. phase', 'Derived surface', 'Derived inter-ply'});
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
        end
        
        function demo_Ascan_plytrack_addnoise(obj, x, y, interp_factor, snr_value)
            % demonstrate the A scan and the ply track
            % including the functions of interpolation,
            % x, y: the x, y index to select the A scan
            % interp_factor: the interpolation multiple
            % snr: snr of the added noise
            Ascan = squeeze(obj.img(x, y, :));
            % interpolation
            xp = 1: 1: length(Ascan);
            xq = 1: 1 / interp_factor: length(Ascan);
            Ascan_interp = interp1(xp, Ascan, xq);
            front_I_interp = obj.front_I(x, y) * interp_factor;
            rear_I_interp = obj.rear_I(x, y) * interp_factor;
            % add noise
            Ascan_interp = awgn(Ascan_interp, snr_value, 'measured');
            % calcualte the snr
            %             SNR = snr(real(Ascan_interp), obj.fs);
            pn = bandpower(Ascan_interp, obj.fs, [10e6, round(obj.fs / 2.1)]);
            ps = bandpower(Ascan_interp, obj.fs, [1e6, 10e6]);
            SNR = 10 * log10(ps / pn);
            disp("SNR:");
            disp(SNR);
            % plot the inam and the record signal
            Ascan_as = hilbert(Ascan_interp);
            cf = figure('Name', strcat('plytrack', num2str(x)));
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ca = subplot(2, 1, 1);
            plot(ca, (1:length(Ascan_interp)) / obj.fs, Ascan_interp, 'linewidth', 2);
            hold on;
            plot(ca, (1:length(Ascan_as)) / obj.fs, abs(Ascan_as), 'linewidth', 2)
            ylabel(ca,{'\fontname {times new roman}\fontsize {16}Amplitude (arb.)'});
            xlabel(ca, '\fontname {times new roman} Time(s)', 'fontsize', 16);
            legend({'Recorded', 'Inst. Amp.'});
            legend('Location', 'best');
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % track the - pi / 2 by interplation on the unwraped data
            A_inph = angle(Ascan_as(round(front_I_interp): round(rear_I_interp))); % inph in the surface window
            A_inph_unw = unwrap(A_inph);
            phase_seq = A_inph_unw(1) + 3 * pi / 2 :  2 * pi : A_inph_unw(end);
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
        
        function show_oneinterply(obj, layer, filtername, v)
            % demonstrate the profile of one inter-ply in 3d
            % layer: 'layer'th layer
            % filtername: 'nofilter' or 'logGabr', for naming the figure.
            % v: the average velocity of the laminate
            % ***
            % calculate the inph
            % The sequence is changed unexpectedly. Transverse the matrix
            inph_ex = obj.mask_plytrack;
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
            ax = subplot(1, 1, 1);
            X = (1: size(profile_layer, 1)) / obj.fx * 1e3;
            Y = (1: size(profile_layer, 2)) / obj.fy * 1e3;
            % convert to distance um from TOF um
            profile_layer = profile_layer * v / 2 / 1e3;
            imagesc(ax, X, Y, profile_layer)
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
            x           = (1:(win(end)-win(1)+1)) / obj.fx * 1e3;
            y           = (1:(win(end)-win(1)+1)) / obj.fx * 1e3;
            z           = (1: size(inph_ex, 3) ) / obj.fs * 1e6;
            inph_visual = abs(obj.(PropertyName));
            % %             mask to select the interply track;
            %             for i = 1:size(inph_ex, 1)
            %                 for j = 1:size(inph_ex, 2)
            %                     indx(i, j, :) = find(inph_ex(i, j, :));
            %                 end
            %             end
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure('Name', ['B_scan_interplytrack_' PropertyName, '_', B_type, '_', num2str(index)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            %             ax = subplot(1, 1, 1);
            if (B_type == 'x')
                % B scan image
                B_scan = squeeze(inph_visual(win, index, :));
                % front line
                front_line = squeeze(obj.front_I(win, index));
                % rear line
                rear_line = squeeze(obj.rear_I(win, index));
                % interply track
                inph_ex = squeeze(inph_ex(win, index, :));
            elseif (B_type == 'y')
                % B scan image
                B_scan = squeeze(inph_visual(index, win, :));
                % front line
                front_line = squeeze(obj.front_I(index, win));
                % rear line
                rear_line = squeeze(obj.rear_I(index, win));
                % interply track
                inph_ex = squeeze(inph_ex(index, win, :));
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
            imagesc(x, z, B_scan');
            hold on;
            % figure setup
            colormap(gray);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16}Inst. amp. (arb)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            % draw the line indicating the front surface
            hold on;
            plot(x, front_line, 'r-', 'linewidth', 1);
            % draw the line indicating the rear surface
            hold on;
            plot(x, rear_line, 'magenta-', 'linewidth', 1);
            hold on;
            scatter(row_inph, col_inph, 2, 'cyano', 'filled');
            legend({'Front surface', 'Rear surface', 'Interply track'});
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
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
        
        function thickness_estimation(obj, total_thick, n)
            % calculate the thickness of each ply
            % demonstrate the thickness by mean and variance
            % total_thick: the total thickness of the composite unit: mm
            % n: the number of the layers to be found
            inph_ex = obj.mask_plytrack;
            thickness = NaN(size(inph_ex, 1), size(inph_ex, 2), n + 1);
            for i = 1: size(inph_ex, 1)
                for j = 1:size(inph_ex, 2)
                    ls =  obj.rear_I(i, j) -  obj.front_I(i, j); % total samples
                    k = find(inph_ex(i, j, :), n);  % if the non-zero values are less that layer, return the real numbers of non-zero values
                    % 1st
                    if k(1) > obj.front_I(i, j)
                        thickness(i, j, 1) = (k(1) - obj.front_I(i, j)) / ls * total_thick ;
                    end
                    thickness(i, j, 2:length(k)) = diff(k) / ls * total_thick ;
                    % last
                    if obj.rear_I(i, j) > k(end)
                        thickness(i, j, n + 1) = (obj.rear_I(i, j) - k(end)) / ls * total_thick ;
                    end
                end
            end
            mean_thickness = squeeze(mean(thickness, [1, 2], 'omitnan'));
            var_thickness = squeeze(var(thickness, 0, [1, 2], 'omitnan'));
            % plot the thickness by mean and var
            figure('Name', ['thickness_estimation_' , num2str(total_thick * 100), '_', num2str(n)]);
            ax = subplot(1, 1, 1);
            errorbar(1:(n + 1), mean_thickness, var_thickness, ...
                '.', 'MarkerSize', 20, 'CapSize', 10, ...
                'MarkerEdgeColor', 'red', 'MarkerFaceColor', 'red');
            hold on;
            grid on;
            xticks(1:2:24);
            xlabel('\fontname {times new roman} layer i^{ply} ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Thickness (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            box on;
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
            % B scan image
            img_WienerDeconv_temp    = NaN(length(B_win), lz);
            img_WienerDeconv_AR_temp = NaN(length(B_win), lz);
            img_SparseDeconv_temp    = NaN(length(B_win), lz);
            if (B_type == 'x')
                % B scan image
                for i = B_win
                    ori_signal = squeeze(obj.img(i, index, :));
                    % remove the direct component
                    ori_signal = ori_signal - mean(ori_signal);
                    % wiener deconvolve
                    wienerdeconv_s                                  = fx_wiener_deconvolution_1d(ori_signal', kernel', q_factor);
                    img_WienerDeconv_temp(i - min(B_win) + 1, :)    = abs(hilbert(wienerdeconv_s));
                    % Wiener deconv. + ARregression
                    deconv = class_deconv(ori_signal, kernel, obj.fs);
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
                    disp([num2str(i - min(B_win) + 1) '/' num2str(lx)]);
                end
            elseif (B_type == 'y')
                for i = B_win
                    ori_signal = squeeze(obj.img(index, i, :));
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
            ca = subplot(3, 1, 2);
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
            wienerdeconv_ARextrap_s     = deconv.wienerdeconv_ARextrap(q_factor_AR, bw, k, 'bg', fft_padding);
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
            caxis(cl / 3);
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
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Magnitude (arb.)\newline ');
            cf = figure('Name', ['Bscan_SparseDeconv' '_', B_type, '_', num2str(index), fig_subname]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            imagesc(x, z, img_SparseDeconv_temp');
            % figure setup
            colormap(gray);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Magnitude (arb.)');
            set(gca, 'Fontname', 'times new Roman');
            xlabel('\fontname {times new roman} X displacement (mm)', 'fontsize', 16);
            ylabel('\fontname {times new roman} TOF (\mus)', 'fontsize', 16);
            set(gca, 'fontsize', 16);
            set(gca, 'linewidth', 2);
            caxis('auto');
            cl =  caxis;
            caxis(cl / 5);
        end
        
        % ************************* In-plane direction *********
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
            inam = abs(obj.(PropertyName));
            inph_ex = obj.mask_plytrack;
            profile_layer = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
            obj.C_scan_inam = NaN([size(inph_ex, 1), size(inph_ex, 2)]);
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
                        obj.C_scan_inam(i, j) = inam(i, j, z_index);
                        profile_layer(i, j) = z_index / obj.fs * 1e6;
                    else % set the last inter-ply as NaN otherwise
                        profile_layer(i, j) = NaN;
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
            Cscan_inam = NaN(size(profile_layer));
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
                        Cscan_inam(i, j) = inam(i, j, z_index);
                        profile_layer_z(i, j) = z_index;
                    end
                end
            end
            % fillna
            Cscan_inam(isnan(Cscan_inam)) = mean(Cscan_inam, 'all', 'omitnan');
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
        
        function obj = define_parallel_inamCscan(obj, dis_to_front, PropertyName, center, radii)
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
            % creat a mask for rear surface
            obj.rear_mask = obj.front_I + dis_to_front < obj.rear_I;
            for i = 1: size(inam, 1)
                for j = 1:size(inam, 2)
                    front_index = obj.front_I(i, j);
                    if front_index + dis_to_front < obj.rear_I(i, j)
                        obj.C_scan_inam(i, j) = inam(i, j, round(front_index + dis_to_front));
                    end
                end
            end
            disp(dis_to_front);
            % creat a round mask
            mask = fx_createCirclesMask(obj.C_scan_inam, center, radii);
            x1 = center(2) - radii;
            x2 = center(2) + radii;
            y1 = center(1) - radii;
            y2 = center(1) + radii;
            % use the mask to select a round region
            C_scan_inam_mask = obj.C_scan_inam .* mask;
            C_scan_inam_mask = C_scan_inam_mask(x1: x2,y1: y2);
            figure('Name', ['Cscan_RT_parallel_' , PropertyName, '_', num2str(dis_to_front)]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_mask, 2)) / obj.fx * 1e3;
            Y = (1: size(C_scan_inam_mask, 1)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_mask);
            hold on;
            h = colorbar;
            % caxis([4.5 5.5]); % for comparision, the 'caxis' needs to be modified
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} inst. amp. (arb,)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            %
            figure('Name', ['Cscan_RT_parallel_' , PropertyName, '_', num2str(dis_to_front)]);
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
                        [~, anguler_1D, ~, ~] = fx_Radonto1Dangular(C_scan_inam_RT_pad, center, r, theta);
                        ID(i, j, r_idx, :) = anguler_1D;
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
        
        function obj = compute_orientation_by_RT_correct(obj, radiis, theta)
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
            C_scan_inam_RT_pad = padarray(C_scan_inam_RT, [r_max r_max], mean(C_scan_inam_RT, 'all', 'omitnan'));
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
                        ID(i, j, r_idx, :) = anguler_1D;
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
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
        end
        
        function show_orientation_by_ID_RT(obj, radius, theta)
            % demonstrate the orientation extracted by Information diagram
            % of RT
            % use the inam image calcualated before
            % radius, theta: the same ratio for the multi-resolution RT
            % ***
            [lx, ly] = size(obj.C_scan_inam);
            % remove the edges
            image_orientation = NaN(lx, ly);
            image_radius = NaN(lx, ly);
            ID_sum = zeros(length(radius), length(theta));
            for i = 1: lx
                for j = 1: ly
                    ID = squeeze(obj.RT_ID(i, j, :, :));
                    [~, I] = max(ID(:));
                    [indrad, indAng] = ind2sub(size(ID),I);
                    image_orientation(i, j) = theta(indAng);
                    image_radius(i, j) = radius(indrad);
                    % sum up the ID
                    ID_sum = ID_sum + ID;
                end
            end
            if ~isempty(obj.rear_mask)
                image_orientation(~obj.rear_mask) = NaN;
            end
            %
            cf = figure();
            ax = subplot(1, 1, 1);
            surf(theta, radius, ID_sum);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Wavelength (pixel)', 'fontsize', 16);
            %
            figure('Name', ['orientation_image_RT_']);
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
        
        function obj = extract_local_orientation_RT_3D_parallel(obj, PropertyName, radius, theta, distances_to_front)
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
                            [~, anguler_1D, ~, ~] = fx_Radonto1Dangular(C_scan_inam_RT_pad, center, r, theta);
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
                        %                         front_index = obj.front_I(ii, jj);
                        %                         % special transverse when one dimension is 1;
                        %                         if size(ID_pixel, 2)==1
                        %                             ID_pixel = ID_pixel';
                        %                         end
                        %                         [~, I] = max(ID_pixel(:));
                        %                         [~, indAng] = ind2sub(size(ID_pixel),I);
                        %                         if front_index + dis_to_front < obj.rear_I(ii, jj) && front_index > 0
                        %                             Inplane_direction(ii, jj, round(front_index + dis_to_front - spacing_dis / 2)...
                        %                                 :round(front_index + dis_to_front + spacing_dis / 2)) = theta(indAng);
                        %                         end
                        % sum up the ID
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
            obj.Inplane_direction_3D = Inplane_direction;
            obj.Inplane_direction_3D_overall = Inplane_direction_overall;
            obj.ID_sum_overall = ID_sum_all;
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
            obj.Inplane_direction_3D_modifyRT = Inplane_direction;
            obj.Inplane_direction_3D_multireso_modifyRT = Inplane_direction_multireso;
            obj.Inplane_direction_3D_overall_modifyRT = Inplane_direction_overall;
            obj.ID_sum_overall_modifyRT = ID_sum_all;
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
                    r = radius(r_idx);
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
        
        function show_inplane_direction_3D(obj, xslice, yslice, zslice)
            % demonstrate the 3d inplane_direction extracted by Gabor filter
            % demonstrate the ID_sum_overall
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % ***
            inph_ex = obj.Inplane_direction_3D;
            % the 1st dimension is length(radius)
            y = (0: size(inph_ex, 2) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 3) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 4) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            for i = 1: size(inph_ex, 1)
                cf = figure('Name', ['3d_orientation_RT_' 'r' num2str(obj.radius_RT(i)) '_' 'xslice', num2str(xslice(1))]);
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
            % ************ Inplane_orientation_3D_overall ******
            inph_ex = obj.Inplane_direction_3D_overall;
            y = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            inph_visual = inph_ex;
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure('Name', ['3d_orientation_RT_overall_' num2str(obj.radius_RT) '_' num2str(xslice(1))]);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            h = slice(ax, X, Y, Z, inph_visual, xslice, yslice, zslice);
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
            ID_sum_all = squeeze(obj.ID_sum_overall);
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
        
        function obj = compute_logGabor_filter_withoutFig(obj, PropertyName, wavelength, orientation)
            % demonstrate and save the logGabor filtered images
            % use the C_scan_inam defined before
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ***
            gaborArray = gabor(wavelength, orientation, 'SpatialAspectRatio', 0.5);
            % fillna
            C_scan_inam_LG = obj.C_scan_inam;
            C_scan_inam_LG(isnan(C_scan_inam_LG)) = mean(C_scan_inam_LG, 'all', 'omitnan');
            % pad the C_scan_inam matrix
            wl_max = max(wavelength);
            %             C_scan_inam_RT_pad = padarray(C_scan_inam_LG, [round(wl_max / 2) round(wl_max / 2)], 'replicate');
            C_scan_inam_RT_pad = padarray(C_scan_inam_LG, [round(wl_max / 2) round(wl_max / 2)], mean(C_scan_inam_LG, 'all', 'omitnan'));
            %
            figure('Name', ['Cscan_fitplytrack_' , PropertyName]);
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan_inam_RT_pad, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan_inam_RT_pad, 2)) / obj.fy * 1e3;
            imagesc(ax, X, Y, C_scan_inam_RT_pad);
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
            % Apply filters to C-scan image. **replication padding**
            [obj.gaborMag, obj.gaborPha] = imgaborfilt(C_scan_inam_RT_pad, gaborArray);
            % normalized the magnitude.
            for i = 1:length(gaborArray)
                BW = gaborArray(i).SpatialFrequencyBandwidth;
                sigmaX = gaborArray(i).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                sigmaY = sigmaX ./ gaborArray(i).SpatialAspectRatio;
                obj.gaborMag(:, :, i) = obj.gaborMag(:, :, i)  / (2 * sigmaX * sigmaY * pi);
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
                    ID = squeeze(obj.gaborMag(i + round(wl_max / 2), j + round(wl_max / 2), :));
                    ID = reshape(ID, [obj.frow, obj.fcol]);
                    angular_1D = sum(ID, 1);
                    % shift the angular_1D by 90 degrees (half length) in x axis
                    angular_1D_shift = circshift(angular_1D, round(length(angular_1D)/2));
                    angular_1D_sub = angular_1D - angular_1D_shift;
                    [~, indAng] = max(angular_1D_sub);
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
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
        end
        
        function show_orientation_by_ID_allwl(obj, wavelength, orientation)
            % demonstrate the orientation extracted by Information diagram
            % use the logGabor filtered images calcualated before
            % orientation: the orientation for Gabor filter, to
            % determine the angle space
            % wavelength: the wavelength for Gabor filter, to determing the size of kernel for removing the edges
            % ***
            % kernel_size_half =
            [lx, ly] = size(obj.C_scan_inam);
            % remove the edges
            wl_max = max(wavelength);
            image_orientation = NaN(lx, ly);
            image_wavelength = NaN(lx, ly);
            ID_sum = zeros(length(wavelength), length(orientation));
            for i = 1: lx
                for j = 1: ly
                    ID = squeeze(obj.gaborMag(i + round(wl_max / 2), j + round(wl_max / 2), :));
                    ID = reshape(ID, [obj.frow, obj.fcol]);
                    %                     ID_shift = circshift(ID, round(size(ID, 2) / 2), 2);
                    %                     ID_sub = ID - ID_shift;
                    [~, I] = max(ID(:));
                    [indwl, indAng] = ind2sub(size(ID),I);
                    image_orientation(i, j) = orientation(indAng);
                    image_wavelength(i, j) = wavelength(indwl);
                    % sum up ID
                    ID_sum = ID_sum + ID;
                end
            end
            %             if ~isempty(obj.rear_mask)
            %                 image_orientation(~obj.rear_mask) = NaN;
            %             end
            %
            cf = figure();
            ax = subplot(1, 1, 1);
            surf(orientation, wavelength, ID_sum);
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Wavelength (pixel)', 'fontsize', 16);
            %
            figure('Name', ['orientation_image']);
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
            %
            figure('Name', ['wavelength_image']);
            ax = subplot(1, 1, 1);
            X = (1: size(image_wavelength, 2)); % pixels
            Y = (1: size(image_wavelength, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
            pcolor(ax, X, Y, image_wavelength);
            set(gca, 'YDir', 'reverse');
            shading flat;
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Wavelength (pixel)');
        end
        
        function show_orientation_by_ID_onewl(obj, wl_idx, wavelength, orientation)
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
        
        function [max_angle_I, obj] = depth_Radontransform_vitualprofile(obj, ratios, PropertyName, center, radii)
            % demonstrate the profile of one inter-ply in 3d
            % save the C_scan_inam
            % layer: number of the layer, 'layer'th layer
            % ratio: the ratio of the position in one ply: 10%, 20%,
            % .....90%
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % ***
            inam = abs(obj.(PropertyName));
            % first layer: front surface
            profile_layer_pre = obj.front_I;
            % loop for depth
            max_angle_I = NaN(1, max(obj.layer_counts_new) * length(ratios));
            m = 1; % count for the max_angle_I
            for layer = 2: max(obj.layer_counts_new)
                % 'layer'th interply
                oneply_row = obj.row(obj.layer_counts_new==layer);
                oneply_col = obj.col(obj.layer_counts_new==layer);  % / obj.fy * 1e3 ;
                oneply_dep = obj.dep(obj.layer_counts_new==layer);  % / obj.fy * 1e3 ;
                %
                profile_layer = NaN(max(oneply_row), max(oneply_col));
                for i = 1: length(oneply_row)
                    profile_layer(oneply_row(i), oneply_col(i)) =  oneply_dep(i);
                end
                % construct the target profile and the inam
                C_scan_inam_depth = NaN(size(profile_layer));
                profile_layer_ave = mean(profile_layer, 'all', 'omitnan');
                profile_layer_pre_ave = mean(profile_layer_pre, 'all', 'omitnan');
                for ratio = 1: length(ratios)
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
                            for ratio = 1: length(ratios)
                                if k2 > k1 && obj.rear_I(i, j) > k2
                                    % layer is deeper than previous layer & rear_surface is deeper than all
                                    z_index = round(k1 + (k2 - k1) * ratio);
                                    C_scan_inam_depth(i, j) = inam(i, j, z_index);
                                end
                            end
                        end
                    end
                end
                C_scan_inam_depth(isnan(obj.C_scan_inam)) = mean(C_scan_inam_depth, 'all', 'omitnan');
                % Radon transform
                % creat a round mask parameters
                [~, ~, max_angle_I(m), ~, ~] = fx_Radonto1Dangular(C_scan_inam_depth, center, radii);
                disp(m);
                m = m + 1;
                profile_layer_pre = profile_layer;
            end
            
        end
        
        function obj = extract_local_orientation_3D(obj, PropertyName, wavelength, orientation, ratios, max_layer)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ratios: the depth in each ply
            % max_layer: the total plies of the structure
            % ***
            gaborArray = gabor(wavelength, orientation);
            row_ga = length(wavelength);
            col_ga = length(orientation);
            Inplane_direction = NaN(size(obj.(PropertyName)));
            ratio_spacing = (ratios(end) - ratios(1)) / (length(ratios) - 1);
            for i = 1:max_layer
                for j = 1:length(ratios)
                    ratio = ratios(j);
                    [Cscan_inam, profile_layer, profile_layer_pre] = obj.define_Cscan_depthprofile_knn(i, ratio, PropertyName, max_layer);
                    % Apply filters to C-scan image.
                    [gaborMagnitude, ~] = imgaborfilt(Cscan_inam, gaborArray);
                    % normalized the magnitude.
                    for k = 1:length(gaborArray)
                        BW = gaborArray(k).SpatialFrequencyBandwidth;
                        sigmaX = gaborArray(k).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                        sigmaY = sigmaX ./ gaborArray(k).SpatialAspectRatio;
                        gaborMagnitude(:, :, k) = gaborMagnitude(:, :, k)  / (2 * sigmaX * sigmaY * pi);
                    end
                    % extract the local orientaion
                    [lx, ly] = size(gaborMagnitude(:, :, 1));
                    image_orientation = zeros(lx, ly);
                    for ii = 1:lx
                        for jj = 1:ly
                            ID = squeeze(gaborMagnitude(ii, jj, :));
                            ID = reshape(ID, [row_ga, col_ga]);
                            angular_1D = sum(ID, 1);
                            % shift the angular_1D by 90 degrees (half length) in x axis
                            angular_1D_shift = circshift(angular_1D, length(angular_1D) / 2);
                            angular_1D_sub = angular_1D - angular_1D_shift;
                            [~, indAng] = max(angular_1D_sub);
                            % save the orientation image
                            if profile_layer_pre(ii, jj) <= profile_layer(ii, jj)
                                h = profile_layer(ii, jj) - profile_layer_pre(ii, jj);
                                k1 = round(profile_layer_pre(ii, jj) + h * (ratio - ratio_spacing));
                                k2 = round(profile_layer_pre(ii, jj) + h * (ratio + ratio_spacing));
                                Inplane_direction(ii, jj, k1:k2) = orientation(indAng);
                                %                                 image_orientation(ii, jj) = orientation(indAng);
                            end
                        end
                    end
                    % end loop in one Cscan
                end
                %                 % debug
                %                 if i==3
                %                     figure('Name', ['orientation_image']);
                %                     ax = subplot(1, 1, 1);
                %                     X = (1: size(image_orientation, 2)); % pixels
                %                     Y = (1: size(image_orientation, 1)); % default X = (1: size(C_scan_inam_LG, 2))  here
                %                     imagesc(ax, X, Y, image_orientation);
                %                     h = colorbar;
                %                     set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Angle (rad.)');
                %                 end
                disp([num2str(i), '/', num2str(max_layer), ' layer']);
            end
            obj.Inplane_direction_3D = Inplane_direction;
        end
        
        function obj = extract_local_orientation_3D_parallel(obj, PropertyName, wavelength, orientation, distances_to_front)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % distances_to_front: the distance to the front surface , unit: us
            % ***
            gaborArray = gabor(wavelength, orientation);
            row_ga = length(wavelength);
            col_ga = length(orientation);
            wl_max = max(wavelength);
            Inplane_direction = NaN(size(obj.(PropertyName)));
            inam = abs(obj.(PropertyName));
            spacing_dis = (distances_to_front(end) - distances_to_front(1)) ...
                / (length(distances_to_front) - 1) * obj.fs / 1e6;
            for i = 1:length(distances_to_front)
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
                % fillna
                C_scan_inam_para(isnan(C_scan_inam_para)) = mean(C_scan_inam_para, 'all', 'omitnan');
                [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_para, gaborArray);
                % normalized the magnitude.
                for k = 1:length(gaborArray)
                    BW = gaborArray(k).SpatialFrequencyBandwidth;
                    sigmaX = gaborArray(k).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                    sigmaY = sigmaX ./ gaborArray(k).SpatialAspectRatio;
                    gaborMagnitude(:, :, k) = gaborMagnitude(:, :, k)  / (2 * sigmaX * sigmaY * pi);
                end
                % extract the local orientaion
                [lx, ly] = size(gaborMagnitude(:, :, 1));
                for ii = round(wl_max / 2) + 1: round(lx - wl_max / 2)
                    for jj = round(wl_max / 2) + 1: round(ly - wl_max / 2)
                        if obj.front_I(ii, jj) + dis_to_front < obj.rear_I(ii, jj)
                            ID = squeeze(gaborMagnitude(ii, jj, :));
                            ID = reshape(ID, [row_ga, col_ga]);
                            [~, I] = max(ID(:));
                            [~, indAng] = ind2sub(size(ID),I);
                            Inplane_direction(ii, jj, round(front_index + dis_to_front - spacing_dis / 2) ...
                                : round(front_index + dis_to_front + spacing_dis / 2)) = orientation(indAng);
                        end
                    end
                end
                % end loop in one Cscan
                disp([num2str(distances_to_front(i)), '/', num2str(distances_to_front(end)), 'us']);
            end
            obj.Inplane_direction_3D = Inplane_direction;
        end
        
        function obj = extract_local_orientation_3D_parallel_allwl(obj, PropertyName, wavelength, orientation, distances_to_front)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % distances_to_front: the distance to the front surface , unit: us
            % ***
            obj.wavelength_LG = wavelength;
            obj.theta_LG = orientation;
            obj.distances_to_front_LG = distances_to_front;
            obj.PropertyName_LG = PropertyName;
            % start
            gaborArray = gabor(wavelength, orientation);
            row_ga = length(wavelength);
            col_ga = length(orientation);
            wl_max = max(wavelength);
            Inplane_direction = NaN(size(obj.(PropertyName)));
            Inplane_direction_overall = NaN(size(obj.(PropertyName)));
            inam = abs(obj.(PropertyName));
            spacing_dis = (distances_to_front(end) - distances_to_front(1)) ...
                / (length(distances_to_front) - 1) * obj.fs / 1e6;
            ID_sum_all = NaN(length(distances_to_front), length(wavelength), length(orientation));
            for i = 1:length(distances_to_front)
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
                % fillna
                C_scan_inam_para(isnan(C_scan_inam_para)) = mean(C_scan_inam_para, 'all', 'omitnan');
                [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_para, gaborArray);
                % normalized the magnitude.
                for k = 1:length(gaborArray)
                    BW = gaborArray(k).SpatialFrequencyBandwidth;
                    sigmaX = gaborArray(k).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                    sigmaY = sigmaX ./ gaborArray(k).SpatialAspectRatio;
                    gaborMagnitude(:, :, k) = gaborMagnitude(:, :, k)  / (2 * sigmaX * sigmaY * pi);
                end
                % extract the local orientaion
                % search for the max
                ID_sum = zeros(length(wavelength), length(orientation));
                [lx, ly] = size(gaborMagnitude(:, :, 1));
                for ii = round(wl_max / 2) + 1: round(lx - wl_max / 2)
                    for jj = round(wl_max / 2) + 1: round(ly - wl_max / 2)
                        if obj.front_I(ii, jj) + dis_to_front < obj.rear_I(ii, jj)
                            ID = squeeze(gaborMagnitude(ii, jj, :));
                            ID = reshape(ID, [row_ga, col_ga]);
                            [~, I] = max(ID(:));
                            [~, indAng] = ind2sub(size(ID),I);
                            Inplane_direction(ii, jj, round(front_index + dis_to_front - spacing_dis / 2) ...
                                : round(front_index + dis_to_front + spacing_dis / 2)) = orientation(indAng);
                        end
                        % sum up the ID
                        ID_sum = ID_sum + ID;
                    end
                end
                % save the summary ID
                ID_sum_all(i, :, :) = ID_sum;
                [~, I] = max(ID_sum(:));
                [~, indAng] = ind2sub(size(ID),I);
                Inplane_direction_overall(:, :, round(front_index + dis_to_front - spacing_dis / 2) ...
                    : round(front_index + dis_to_front + spacing_dis / 2)) = orientation(indAng);
                % end loop in one Cscan
                disp([num2str(distances_to_front(i)), '/', num2str(distances_to_front(end)), 'us']);
            end
            obj.Inplane_direction_3D_ID = Inplane_direction;
            obj.Inplane_direction_3D_overall_ID = Inplane_direction_overall;
            obj.ID_sum_overall_ID = ID_sum_all;
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
                        [fileID, ~] = fopen(fullFileName, 'w');
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
                        ID_d = squeeze(gaborMag_d(ii + round(wl_max), jj + round(wl_max), :));
                        ID_shape = reshape(ID_d, [lgrow, lgcol]);
                        baseFileName = ['x_' num2str(ii), '_', 'y_', num2str(jj), '_ID.txt'];
                        fullFileName = fullfile(folder, baseFileName);
                        [fileID, ~] = fopen(fullFileName, 'a');
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
        
        % not correct yet
        function obj = extract_local_orientation_3D_byuwp_allwl(obj, PropertyName, wavelength, orientation, distances_to_front)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % distances_to_front: the distance to the front surface , unit: us
            % ***
            gaborArray = gabor(wavelength, orientation);
            row_ga = length(wavelength);
            col_ga = length(orientation);
            wl_max = max(wavelength);
            Inplane_direction = NaN(size(obj.(PropertyName)));
            spacing_dis = (distances_to_front(end) - distances_to_front(1)) ...
                / (length(distances_to_front) - 1) * obj.fs / 1e6;
            for i = 1:length(distances_to_front)
                dis_to_front = distances_to_front(i) * obj.fs / 1e6; % us to sample points
                [Cscan_inam_uwp, obj] = obj.define_Cscan_depthprofile_uwphase(dis_to_front, PropertyName);
                % fillna
                Cscan_inam_uwp(isnan(Cscan_inam_uwp)) = mean(Cscan_inam_uwp, 'all', 'omitnan');
                %                 [gaborMagnitude, ~] = imgaborfilt(Cscan_inam_uwp, gaborArray);
                %                 % normalized the magnitude.
                %                 for k = 1:length(gaborArray)
                %                     BW = gaborArray(k).SpatialFrequencyBandwidth;
                %                     sigmaX = gaborArray(k).Wavelength / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
                %                     sigmaY = sigmaX ./ gaborArray(k).SpatialAspectRatio;
                %                     gaborMagnitude(:, :, k) = gaborMagnitude(:, :, k)  / (2 * sigmaX * sigmaY * pi);
                %                 end
                %                 % extract the local orientaion
                %                 [lx, ly] = size(gaborMagnitude(:, :, 1));
                %                 for ii = round(wl_max / 2) + 1: round(lx - wl_max / 2)
                %                     for jj = round(wl_max / 2) + 1: round(ly - wl_max / 2)
                %                         if obj.front_I(ii, jj) + dis_to_front < obj.rear_I(ii, jj)
                %                             ID = squeeze(gaborMagnitude(ii, jj, :));
                %                             ID = reshape(ID, [row_ga, col_ga]);
                %                             [~, I] = max(ID(:));
                %                             [~, indAng] = ind2sub(size(ID),I);
                %                             Inplane_direction(ii, jj, round(front_index + dis_to_front - spacing_dis / 2) ...
                %                                 : round(front_index + dis_to_front + spacing_dis / 2)) = orientation(indAng);
                %                         end
                %                     end
                %                 end
                % end loop in one Cscan
                disp([num2str(distances_to_front(i)), '/', num2str(distances_to_front(end)), 'us']);
            end
            obj.Inplane_direction_3D = Inplane_direction;
        end
        
        function show_inplane_direction_3D_LG(obj, xslice, yslice, zslice)
            % demonstrate the 3d inplane_direction extracted by Gabor filter
            % demonstrate the ID_sum_overall
            % xslice, yslice, zslice: the index to select the slices,
            % units: mm, mm, us
            % ***
            inph_ex = obj.Inplane_direction_3D_ID;
            y = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure('Name', ['3d_orientation_LG_' num2str(obj.wavelength_LG) '_xslice' num2str(xslice(1))]);
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
            inph_ex = obj.Inplane_direction_3D_overall_ID;
            y = (0: size(inph_ex, 1) - 1) / obj.fx * 1e3;
            x = (0: size(inph_ex, 2) - 1) / obj.fy * 1e3;
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure('Name', ...
                ['3d_orientation_LG_overall' num2str(obj.wavelength_LG(1)) 'to' ...
                num2str(obj.wavelength_LG(end)) num2str(xslice(1))]);
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
            ID_sum_all = squeeze(obj.ID_sum_overall_ID);
            for i = 1:size(ID_sum_all, 1)
                ID_sum_all(i, :) = ID_sum_all(i, :) / max(ID_sum_all(i, :));
            end
            cf = figure('Name', ['3d_ID_sum_overall' '_']);
            set(cf, 'Position', [0, 0, 800, 600], 'color', 'white');
            ax = subplot(1, 1, 1);
            imagesc(ax, squeeze(sum(ID_sum_all, 2)));
            hold on;
            set(h, 'EdgeColor', 'none');
            colormap jet;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Mag. of angular distribution from ID (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} Angle (degree) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Depth (pixel)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
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
            z = (0: size(inph_ex, 3) - 1) / obj.fs * 1e6;
            [X, Y, Z] = meshgrid(x, y, z);
            % inph_visual(mask_phase_interp_visual==1) = NaN;
            cf = figure('Name', ['3dtrack_interply_' num2str(xslice(1))]);
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
            zlabel('\fontname {times new roman} Time(\mus)', 'fontsize', 16);
            set(gca, 'ZDir', 'reverse');
            % select the xslice
            x_idx = xslice * obj.fx / 1e3;
            y_idx = yslice * obj.fy / 1e3;
            mask = (obj.col==x_idx) | (obj.row==y_idx);
            x_row = obj.row(mask) / obj.fx * 1e3;
            x_col = obj.col(mask) / obj.fy * 1e3;  % / obj.fy * 1e3 ;
            x_dep = obj.dep(mask) / obj.fs * 1e6;  % / obj.fy * 1e3 ;
            scatter3(ax, x_col, x_row, x_dep, ...
                3, 'cyan', 'filled', ...
                'DisplayName','Interply track');
            hold on;
            scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                y, obj.front_I(:, x_idx)/ obj.fs * 1e6, ...
                3, 'red', 'filled', ...
                'DisplayName','Front surface');
            hold on;
            scatter3(xslice * ones(1, size(inph_ex, 1)), ...
                y, obj.rear_I(:, x_idx)/ obj.fs * 1e6, ...
                3, 'magenta', 'filled', ...
                'DisplayName','Rear surface');
            hold on;
            % select the yslice
            scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                obj.front_I(y_idx, :)/ obj.fs * 1e6, ...
                3, 'red', 'HandleVisibility','off');
            hold on;
            scatter3(x, yslice * ones(1, size(inph_ex, 2)), ...
                obj.rear_I(y_idx, :)/ obj.fs * 1e6, ...
                3, 'magenta', 'HandleVisibility','off');
            hold off;
            legend;
        end
        
        % finish tne methods
    end
    
end


