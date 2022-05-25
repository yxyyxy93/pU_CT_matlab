classdef class_process_woven_RFdata < class_process_RFdata
    % sub class of class_process_RFdata
    properties
        img_3dfft
        img_fftfilter
        % fringerprint of woven sample
        fft2d_abs_fingerprint
        fft2d_pha_fingerprint
        fft2d_mask_fingerprint
        % 3D ID 
        Inplane_direction_3D_ID_ztheta
    end
    
    methods
        function obj = class_process_woven_RFdata(filename)
            %
            obj = obj@class_process_RFdata(filename);
        end
        
        %************************* A_scan analysis *************************
        function signal = demo_AS_3D_inclinfq(obj, x, y)
            % demonstrate analytic-signal
            % x, y: index of the A scan
            ascan   = squeeze(obj.img(x, y, :));
            S_hil   = hilbert(ascan);
            S_real  = real(S_hil);
            S_imag  = imag(S_hil);
            t_space = (1:length(ascan)) / obj.fs;
            % visual for debug
            % decomposed analytical signal
            cf                  = figure('Name', 'Ascan_Analytic_Signal');
            set(cf, 'Position', [0, 0, 800, 1000], 'color', 'white');
            ca = subplot(5, 1, 1);
            plot(ca, t_space, S_real, 'k-', 'Linewidth', 2);
            %                 yl                             = ylim;
            hold on;
            label_h             = xlabel({'\fontname {times new roman}Time (\mus)', '','(a)'}, 'fontsize', 12);
            label_h.Position(2) = label_h.Position(2) + 0.2;
            ylabel(ca, {'\fontname {times new roman} Amplitude (arb.)'}, 'fontsize', 12);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 12);
            set(ca, 'linewidth', 2) ;
            set(ca, 'visible', 'on');
            % adjust pos.
            p       = get(ca, 'Position');
            p(2)    = p(2) + 0.05;
            set(ca, 'Position', p);
            xlim([t_space(1) t_space(end)]);
            %                 title(ca, 'Original signal', 'fontsize', 16, 'fontname', 'times new roman');
            %
            ca = subplot(5, 1, 2);
            plot3(ca, t_space, S_real, S_imag, 'c-', 'Linewidth', 2);
            hold on;
            label_h              = xlabel({'\fontname {times new roman} Time (\mus)'}, 'fontsize', 12);
            label_h.Position(2)  = label_h.Position(2) + 0.2;
            set(get(ca,'xlabel'),'rotation', 2);
            label_h              = ylabel({'\fontname {times new roman} Real amp. (arb.)' }, 'fontsize', 12);
            set(get(ca,'ylabel'),'rotation', -15);
            label_h.Position(1)  = -1; % change horizontal position of ylabel
            label_h.Position(2)  = 2; % change vertical position of ylabel
            zlabel({'\fontname {times new roman} Imag. amp. (arb.)' }, 'fontsize', 12);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 12);
            set(ca, 'linewidth', 2) ;
            set(ca, 'visible', 'on');
            view(ca, [10 -75 -15]);
            % adjust pos.
            p          = get(ca, 'Position');
            p(2)       = p(2) + 0.06;
            set(ca, 'Position', p);
            % phasor representation
            index_plus = round(0.38*obj.fs/1e6); % set manually
            x          = [t_space(index_plus), t_space(index_plus)];
            y          = [0 S_real(index_plus)]; % set manually
            z          = [0 S_imag(index_plus)];
            line(x, y, z, 'Color', 'Red', 'LineStyle', '-',  'Linewidth', 2);
            scatter3(x(end), y(end), z(end), 25,...
                'MarkerEdgeColor', 'red', ...
                'MarkerFaceColor', 'red',...
                'LineWidth', 1.5);
            xlim([t_space(1) t_space(end)]);
            %
            ca = subplot(5, 1, 3);
            plot(ca, t_space, abs(S_hil), '-', 'Linewidth', 2, 'color', obj.inst_amp_color);
            hold on;
            %                 xlabel({'\fontname {times new roman}Time (\mus)', '', '(c)'}, 'fontsize', 12);
            ylabel({'\fontname {times new roman} Amplitude (arb.)' }, 'fontsize', 12);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 12);
            set(ca, 'linewidth', 2) ;
            set(ca, 'visible', 'on');
            xlim([t_space(1) t_space(end)]);
            % adjust pos.
            p          = get(ca, 'Position');
            p(2)       = p(2) + 0.05;
            set(ca, 'Position', p);
            %                 title(ca, 'Instantaneous amplitude', 'fontsize', 16, 'fontname', 'times new roman');
            xticks([]);
            %
            ca = subplot(5, 1, 4);
            plot(ca, t_space, angle(S_hil), '-', 'Linewidth', 2, 'color', obj.inst_phase_color);
            hold on;
            %                 xlabel({'\fontname {times new roman}Time (\mus)', '(d)'}, 'fontsize', 12);
            ylabel({'\fontname {times new roman}Phase (rad.)' }, 'fontsize', 12);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 12);
            set(ca, 'linewidth', 2) ;
            set(ca, 'visible', 'on');
            xlim([t_space(1) t_space(end)]);
            % adjust pos.
            p          = get(ca, 'Position');
            p(2)       = p(2) + 0.05;
            set(ca, 'Position', p);
            xticks([]);
            %
            ca = subplot(5, 1, 5);
            [ifq, t] = instfreq(S_hil, obj.fs, 'Method', 'hilbert');
            plot(ca, t*1e6, ifq/1e6, '-', 'Linewidth', 2, 'color', [0.9100 0.4100 0.1700]);
            hold on;
            xlabel({'\fontname {times new roman}Time (\mus)'}, 'fontsize', 12);
            ylabel({'\fontname {times new roman}Freq. (MHz)' }, 'fontsize', 12);
            set(ca, 'Fontname', 'times new Roman', 'Fontsize', 12);
            set(ca, 'linewidth', 2) ;
            set(ca, 'visible', 'on');
            xlim([t(1) t(end)]*1e6);
            ylim([0 10]);
            % adjust pos.
            p          = get(ca, 'Position');
            p(2)       = p(2) + 0.05;
            set(ca, 'Position', p);
            % output
            signal = S_hil;
        end
        
        function obj = removeslope(obj, PropertyName)
            % remove
            img_temp = obj.(PropertyName);
            img_hilbert = nan(size(img_temp));
            for i = 1:size(img_temp, 1)
                for j = 1:size(img_temp, 2)
                    img_temp(i, j, :) = img_temp(i, j, :) - mean(img_temp(i, j, :));
                    img_hilbert(i, j, :) = hilbert(img_temp(i, j, :));
                end
                clc;
                disp([num2str(i) '/' num2str(size(img_temp, 1))]);
            end
            obj.(PropertyName) = img_temp;
            obj.img_hil = img_hilbert;
        end
        
        %************************* B_scan analysis *************************
        function demo_Bscan_inst(obj, B_type, index, Bwin, PropertyName)
            % demonstrate the B_scan of the img
            % show inph_ex
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % Bwin: the indexes of the first and the last points in Bscan
            % PropertyName: the dataset to use
            % image
            inph_visual = obj.(PropertyName);
            x = (1: Bwin(end)-Bwin(1)+1)  / obj.fx * 1e3;
            y = (1: Bwin(end)-Bwin(1)+1)  / obj.fy * 1e3;
            z = (1: size(inph_visual, 3)) / obj.fs * 1e6;
            %             ax = subplot(1, 1, 1);
            if (B_type == 'x')
                % B scan image
                B_scan = squeeze(inph_visual(Bwin, index, :));
                % surfaecs xslice
                if ~isnan(obj.front)
                    sur_f   = obj.front_I(Bwin, index)/ obj.fs * 1e6;
                    sur_r   = obj.rear_I(Bwin, index)/ obj.fs * 1e6;
                end
                x_slice = index * ones(1, size(inph_visual, 1));
                y_index = x;
            elseif (B_type == 'y')
                % B scan image
                B_scan = squeeze(inph_visual(index, Bwin, :));
                % surfaecs yslice
                if ~isnan(obj.front)
                    sur_f   = obj.front_I(index, Bwin)/ obj.fs * 1e6;
                    sur_r   = obj.rear_I(index, Bwin)/ obj.fs * 1e6;
                end
                x_slice = index * ones(1, size(inph_visual, 2));
                y_index = x;
                % choose the 'x' axis
                x       = y;
            end
            %             %remove the direct component
            %             B_scan = B_scan - mean(B_scan, 2);
            AS_bscan       = B_scan';
            inph_ex        = angle(AS_bscan);
            inam_ex        = abs(AS_bscan);
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
            cf = figure('Name', ['Bscan' '_', B_type, '_', num2str(index), '_inph']);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            imagesc(x, z(2:end), inph_ex);
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
            imagesc(x, z, infq_ex/1e6);
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
        
        function demo_Bscan_inst_needHilbert(obj, B_type, index, Bwin, PropertyName)
            % demonstrate the B_scan of the img
            % show inph_ex
            % index: the index of the slics in 'B_type' direction, units: mm, mm, us
            % Bwin: the indexes of the first and the last points in Bscan
            % PropertyName: the dataset to use
            % image
            inph_visual = obj.(PropertyName);
            x = (1: Bwin(end)-Bwin(1)+1)  / obj.fx * 1e3;
            y = (1: Bwin(end)-Bwin(1)+1)  / obj.fy * 1e3;
            z = (1: size(inph_visual, 3)) / obj.fs * 1e6;
            %             ax = subplot(1, 1, 1);
            if (B_type == 'x')
                % B scan image
                B_scan = squeeze(inph_visual(Bwin, index, :));
                % surfaecs xslice
                if ~isnan(obj.front)
                    sur_f   = obj.front_I(Bwin, index)/ obj.fs * 1e6;
                    sur_r   = obj.rear_I(Bwin, index)/ obj.fs * 1e6;
                end
                x_slice = index * ones(1, size(inph_visual, 1));
                y_index = x;
            elseif (B_type == 'y')
                % B scan image
                B_scan = squeeze(inph_visual(index, Bwin, :));
                % surfaecs yslice
                if ~isnan(obj.front)
                    sur_f   = obj.front_I(index, Bwin)/ obj.fs * 1e6;
                    sur_r   = obj.rear_I(index, Bwin)/ obj.fs * 1e6;
                end
                x_slice = index * ones(1, size(inph_visual, 2));
                y_index = x;
                % choose the 'x' axis
                x       = y;
            end
            % apply hilbert
            for i = 1:size(B_scan, 1)
                B_scan(i, :) = hilbert(B_scan(i, :));
            end
            %             %remove the direct component
            %             B_scan = B_scan - mean(B_scan, 2);
            AS_bscan       = B_scan';
            inph_ex        = angle(AS_bscan);
            inam_ex        = abs(AS_bscan);
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
            cf = figure('Name', ['Bscan' '_', B_type, '_', num2str(index), '_inph']);
            set(cf, 'Position', [0, 0, 600, 600], 'color', 'white');
            imagesc(x, z(2:end), inph_ex);
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
            imagesc(x, z, infq_ex/1e6);
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
        
        %************************* C-scan analysis ************************
        function show_Cscan_2dfft(obj, z, PropertyName)
            % show the C scan by z index in the depth, and the 2d fft
            % z: z index. unit: data points
            % PropertyName: 'img_hil' or 'img_hil_filter' ...
            % In the image, the unit is transfered to us.
            figure('Name', ['C_scan_amp' , '_', num2str(z)]);
            C_scan = abs(obj.(PropertyName)(:, :, z));
            ax = subplot(1, 1, 1);
            X = (1: size(C_scan, 1)) / obj.fx * 1e3;
            Y = (1: size(C_scan, 2)) / obj.fy * 1e3;
            imagesc(ax, Y, X, C_scan / obj.fx * 1e6)
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            set(ax, 'fontsize', 16);
            set(ax, 'linewidth', 1.5);
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            % 2d fft
            DP = fftshift(fft2(C_scan));
            figure,
            imagesc(abs(log10(DP)));
            axis image;
        end
        
        function apply2dfft_LPandHPfilter(obj, z, PropertyName, R_lo_H, R_Lo, roi)
            % apply 2d fft filter (LP & HP) to the C scan by z index.
            % depth.
            % z: z index. unit: data points
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % R_lo_H: the filter size parameter for high pass fitler (1 - Lo_h)
            % R_Lo: the filter size parameter for high pass fitler (1 - Lo_h)
            % roi: x and y index of ROI
            % ***
            % fillna
            C_scan_inam = abs(obj.(PropertyName)(:, :, z));
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            % display snr
            disp('snr and cnr image_origin:')
            [snr, cnr] = fx_image_snr_2(C_scan_inam, roi);
            disp([snr cnr]);
            % ** low pass fitler in frequency domain
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
            [Ymesh, Xmesh] = meshgrid(Y, X);
            %             Lo_h           = exp(-((Xmesh - Cx).^2 + (Ymesh - Cy).^2)./(2 * R_lo_H).^2);
            %             Hi             = 1 - Lo_h; % High pass filter=1-low pass filter
            Lo             = exp(-(Xmesh.^2 + Ymesh.^2)./(2 * R_Lo).^2);
            % **** plot original
            cf = figure('Name', ['C_scan_amp' , '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 2);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            imagesc(ax, Yd, Xd, C_scan_inam);
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
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
            % ** plot filters
            %             figure('Name', ['1-HPfilter' , '_', num2str(z)]);
            %             imagesc(X, Y, Lo_h);
            %             colorbar;
            cf = figure('Name', ['LPfilter' , '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 1);
            imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, 20*log10(abs(fft2D_shifted)/nfftx/nffty),  'AlphaData', Lo); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %2d fft and filtering
            %             fft2D_shifted_filter = fft2D_shifted .* Hi .* Lo;
            fft2D_shifted_filter = fft2D_shifted .* Lo;
            %             fft2D_shifted_filter = fft2D_shifted;
            fft2D_ishift         = ifftshift(fft2D_shifted_filter);
            C_scan_inam_2DFFT    = abs(ifft2(fft2D_ishift));
            % plot 2d fft filtered
            ax     = subplot(2, 1, 2);
            imagesc(Yd, Xd, abs(C_scan_inam_2DFFT(1: m, 1:n)));
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %             imagesc(X, Y, Lo);
            %             colorbar;
            %             figure('Name', ['2dFreqSpectrum' , '_', num2str(z)]);
            %             imagesc(obj.fx*(-nfftx/2+1: nfftx/2)/nfftx, obj.fy*(-nffty/2+1:nffty/2)/nffty, 20*log10(abs(fft2D_shifted)/nfftx/nffty)); % divide by the x and y length: nfftx and nffty
            %             h = colorbar; colormap(jet);
            %             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            %             xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            %             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(gca, 'linewidth', 2);
            %             %2d fft and filtering
            %             %             fft2D_shifted_filter = fft2D_shifted .* Hi .* Lo;
            %             fft2D_shifted_filter = fft2D_shifted .* Lo;
            %             %             fft2D_shifted_filter = fft2D_shifted;
            %             fft2D_ishift         = ifftshift(fft2D_shifted_filter);
            %             C_scan_inam_2DFFT    = ifft2(fft2D_ishift);
            %             % plot 2d fft filtered
            %             figure('Name', ['C_scan_amp_2dfftFiltered' , '_', num2str(z)]);
            %             imagesc(Yd, Xd, C_scan_inam_2DFFT(1: m, 1:n));
            %             axis image;
            %             hold on;
            %             h = colorbar; colormap(jet);
            %             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            %             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            %             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(gca, 'linewidth', 2);
            % display snr
            disp('LPfilter snr and cnr image_filter:');
            [snr, cnr] = fx_image_snr_2(C_scan_inam_2DFFT(1: m, 1:n), roi);
            disp([snr cnr]);
        end
        
        function apply2dfft_HPfilter(obj, z, PropertyName,R_lo_H, roi)
            % apply 2d fft filter (HP) to the C scan by z index.
            % depth.
            % z: z index. unit: data points
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % R_lo_H: the filter size parameter for high pass fitler (1 - Lo_h)
            % roi: x and y index of ROI
            % ***
            % fillna
            C_scan_inam = abs(obj.(PropertyName)(:, :, z));
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            % display snr
            disp('snr and cnr image_origin:')
            [snr, cnr] = fx_image_snr_2(C_scan_inam, roi);
            disp([snr cnr]);
            % ** low pass fitler in frequency domain
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
            [Ymesh, Xmesh] = meshgrid(Y, X);
            Lo_h           = exp(-(Xmesh.^2 + Ymesh.^2)./(2 * R_lo_H).^2);
            Hi             = 1 - Lo_h; % High pass filter=1-low pass filter
            % **** plot original
            cf = figure('Name', ['C_scan_amp' , '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 2);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            imagesc(ax, Yd, Xd, C_scan_inam);
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
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
            % ** plot filters
            cf = figure('Name', ['LPfilter' , '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 1);
            imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, 20*log10(abs(fft2D_shifted)/nfftx/nffty),  'AlphaData', Hi); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %2d fft and filtering
            fft2D_shifted_filter = fft2D_shifted .* Hi;
            fft2D_ishift         = ifftshift(fft2D_shifted_filter);
            C_scan_inam_2DFFT    = abs(ifft2(fft2D_ishift));
            % plot 2d fft filtered
            ax     = subplot(2, 1, 2);
            imagesc(Yd, Xd, abs(C_scan_inam_2DFFT(1: m, 1:n)));
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % display snr
            disp('HPfilter snr and cnr image_filter:');
            [snr, cnr] = fx_image_snr_2(C_scan_inam_2DFFT(1: m, 1:n), roi);
            disp([snr cnr]);
        end
        
        function apply2dfft_LPsinc(obj, z, PropertyName, freq, roi)
            % apply 2d fft filter (LP & HP) to the C scan by z index.
            % depth.
            % z: z index. unit: data points
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % R_lo_H: the filter size parameter for high pass fitler (1 - Lo_h)
            % R_Lo: the filter size parameter for high pass fitler (1 - Lo_h)
            % roi: x and y index of ROI
            % ***
            % fillna
            C_scan_inam = abs(obj.(PropertyName)(:, :, z));
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            % display snr
            disp('snr and cnr image_origin:')
            [snr, cnr] = fx_image_snr_2(C_scan_inam, roi);
            disp([snr cnr]);
            % ** low pass fitler in frequency domain
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
            [Xmesh, Ymesh] = meshgrid(Y, X);
            R_center       = sqrt(Xmesh.^2 + Ymesh.^2);
            % Create 2d sinc function
            Lo               = sin(pi*R_center*freq)./(pi*R_center*freq);
            Lo(end/2, end/2) = 1;
            % **** plot original
            cf = figure('Name', ['C_scan_amp' , '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 2);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            imagesc(ax, Yd, Xd, C_scan_inam);
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
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
            % ** plot filters
            cf = figure('Name', ['LPfilter' , '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 1);
            imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, 20*log10(abs(fft2D_shifted)/nfftx/nffty),  'AlphaData', Lo); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %2d fft and filtering
            %             fft2D_shifted_filter = fft2D_shifted .* Hi .* Lo;
            fft2D_shifted_filter = fft2D_shifted .* Lo;
            %             fft2D_shifted_filter = fft2D_shifted;
            fft2D_ishift         = ifftshift(fft2D_shifted_filter);
            C_scan_inam_2DFFT    = abs(ifft2(fft2D_ishift));
            % plot 2d fft filtered
            ax     = subplot(2, 1, 2);
            imagesc(Yd, Xd, abs(C_scan_inam_2DFFT(1: m, 1:n)));
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % display snr
            disp('LPfilter snr and cnr image_filter:');
            [snr, cnr] = fx_image_snr_2(C_scan_inam_2DFFT(1: m, 1:n), roi);
            disp([snr cnr]);
        end
        
        function apply2dfft_bandreject(obj, z, PropertyName, D_0, W, n_bwf, F_type, roi)
            % apply 2d fft filter (bandreject) to the C scan by z index.
            % depth.
            % z: z index. unit: data points
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % W: width of the band
            % D_0: radius of the circle-shaped band.
            % n_bwf: order of the Butterworth bandreject filter
            % F_type: choose the type of the filter: 'ideal', 'Butterworth', or 'Gaussian'
            % roi: x and y index of ROI
            % ***
            % fillna
            C_scan_inam = abs(obj.(PropertyName)(:, :, z));
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            % display snr
            disp('*********************************************');
            disp('snr and cnr image_origin:')
            [snr, cnr] = fx_image_snr_2(C_scan_inam, roi);
            disp([snr cnr]);
            % ** band reject fitler in frequency domain
            % now make 2D fft of original image
            nfftx         = 2^nextpow2(m);
            nffty         = 2^nextpow2(n);
            % cannot solve the problem when nfftx!=nffty !
            nfftx         = max(nfftx, nffty);
            nffty         = max(nfftx, nffty);
            %
            fft2D         = fft2(C_scan_inam, nfftx, nffty);
            fft2D_shifted = fftshift(fft2D);
            X              = -nfftx/2+1: nfftx/2;
            Y              = -nffty/2+1:nffty/2;
            [Xmesh, Ymesh] = meshgrid(Y, X);
            % design the filter
            D              = sqrt(Xmesh.^2 + Ymesh.^2);
            switch F_type
                case char('ideal')
                    filter         = ones(nfftx, nffty);
                    filter(D<(D_0+W/2) & D>(D_0-W/2)) = 0;
                case char('Butterworth')
                    % Butterworth bandreject filter
                    filter         = 1 ./ (1 + (D*W./(D.^2-D_0^2)).^(2*n_bwf));
                case char('Gaussian')
                    % Gaussian bandreject filter
                    filter         = 1 - exp(-1/2 * ((D.^2-D_0^2)./D./W).^2);
            end
            %
            cf = figure('Name', ['C_scan_amp_2dfftFiltered', '_', F_type, '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 1);
            imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, 20*log10(abs(fft2D_shifted)/nfftx/nffty),  'AlphaData', filter); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %2d fft and filtering
            %             fft2D_shifted_filter = fft2D_shifted .* Hi .* Lo;
            fft2D_shifted_filter = fft2D_shifted .* filter;
            %             fft2D_shifted_filter = fft2D_shifted;
            fft2D_ishift         = ifftshift(fft2D_shifted_filter);
            C_scan_inam_2DFFT    = real(ifft2(fft2D_ishift));
            % plot 2d fft filtered
            ax     = subplot(2, 1, 2);
            imagesc(Yd, Xd, abs(C_scan_inam_2DFFT(1: m, 1:n)));
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % display snr
            disp([F_type 'snr and cnr image_filter:']);
            [snr, cnr] = fx_image_snr_2(C_scan_inam_2DFFT(1: m, 1:n), roi);
            disp([snr cnr]);
        end
        
        function apply2dfft_bandpass(obj, z, PropertyName, D_0, W, n_bwf, F_type, roi)
            % apply 2d fft filter (band pass) to the C scan by z index.
            % depth.
            % z: z index. unit: data points
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % W: width of the band
            % D_0: radius of the circle-shaped band.
            % n_bwf: order of the Butterworth bandreject filter
            % F_type: choose the type of the filter: 'ideal', 'Butterworth', or 'Gaussian'
            % roi: x and y index of ROI
            % ***
            % fillna
            C_scan_inam = abs(obj.(PropertyName)(:, :, z));
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            % display snr
            disp('*********************************************');
            disp('snr and cnr image_origin:')
            [snr, cnr] = fx_image_snr_2(C_scan_inam, roi);
            disp([snr cnr]);
            % ** band reject fitler in frequency domain
            % now make 2D fft of original image
            nfftx         = 2^nextpow2(m);
            nffty         = 2^nextpow2(n);
            % cannot solve the problem when nfftx!=nffty !
            nfftx         = max(nfftx, nffty);
            nffty         = max(nfftx, nffty);
            %
            fft2D         = fft2(C_scan_inam, nfftx, nffty);
            fft2D_shifted = fftshift(fft2D);
            X              = -nfftx/2+1: nfftx/2;
            Y              = -nffty/2+1:nffty/2;
            [Xmesh, Ymesh] = meshgrid(Y, X);
            % design the filter
            D              = sqrt(Xmesh.^2 + Ymesh.^2);
            switch F_type
                case char('ideal')
                    filter         = ones(nfftx, nffty);
                    filter(D>(D_0+W/2) | D<(D_0-W/2)) = 0;
                case char('Butterworth')
                    % Butterworth bandreject filter
                    filter         = 1 - 1 ./ (1 + (D*W./(D.^2-D_0^2)).^(2*n_bwf));
                case char('Gaussian')
                    % Gaussian bandreject filter
                    filter         = exp(-1/2 * ((D.^2-D_0^2)./D./W).^2);
            end
            %
            cf = figure('Name', ['C_scan_amp_2dfftFiltered', '_', F_type, '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 1);
            imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, 20*log10(abs(fft2D_shifted)/nfftx/nffty),  'AlphaData', filter); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %2d fft and filtering
            %             fft2D_shifted_filter = fft2D_shifted .* Hi .* Lo;
            fft2D_shifted_filter = fft2D_shifted .* filter;
            %             fft2D_shifted_filter = fft2D_shifted;
            fft2D_ishift         = ifftshift(fft2D_shifted_filter);
            C_scan_inam_2DFFT    = real(ifft2(fft2D_ishift));
            % plot 2d fft filtered
            ax     = subplot(2, 1, 2);
            imagesc(Yd, Xd, abs(C_scan_inam_2DFFT(1: m, 1:n)));
            axis image;
            hold on;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % display snr
            disp([F_type 'snr and cnr image_filter:']);
            [snr, cnr] = fx_image_snr_2(C_scan_inam_2DFFT(1: m, 1:n), roi);
            disp([snr cnr]);
        end
        
        function apply2dfft_notch(obj, z, PropertyName, D_0, Xv_0, Yv_0, n_bwf, F_type, roi)
            % apply 2d fft filter (notch) to the C scan by z index.
            % z: z index. unit: data points
            % center: [j, i], the centers of the mask for Radon transform
            % radii: the radii of of the mask for Radon transform
            % X_0, Y_0: centers (X_0, Y_0) of the notches
            % D_0: radius of the notches.
            % n_bwf: order of the Butterworth bandreject filter
            % F_type: choose the type of the filter: 'ideal', 'Butterworth', or 'Gaussian'
            % roi: x and y index of ROI
            % ***
            % fillna
            C_scan_inam = abs(obj.(PropertyName)(:, :, z));
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            % plot original
            % display snr
            disp('*************************************************')
            disp('snr and cnr image_origin:')
            [snr, cnr] = fx_image_snr_2(C_scan_inam, roi);
            disp([snr cnr]);
            % ** low pass fitler in frequency domain
            % now make 2D fft of original image
            nfftx         = 2^nextpow2(m);
            nffty         = 2^nextpow2(n);
            % cannot solve the problem when nfftx!=nffty !
            nfftx         = max(nfftx, nffty);
            nffty         = max(nfftx, nffty);
            %
            fft2D         = fft2(C_scan_inam, nfftx, nffty);
            fft2D_shifted = fftshift(fft2D);
            X              = -nfftx/2+1: nfftx/2;
            Y              = -nffty/2+1:nffty/2;
            [Xmesh, Ymesh] = meshgrid(Y, X);
            % design the filter
            fft2D_shifted_filter = fft2D_shifted;
            for i = 1:length(Xv_0)
                X_0 = Xv_0(i);
                Y_0 = Yv_0(i);
                D1  = sqrt((Xmesh-X_0).^2 + (Ymesh-Y_0).^2);
                D2  = sqrt((Xmesh+X_0).^2 + (Ymesh+Y_0).^2); % symmetry
                switch F_type
                    case char('ideal')
                        filter         = ones(nfftx, nffty);
                        filter(D1<D_0 | D2<D_0) = 0;
                    case char('Butterworth')
                        % Butterworth bandreject filter
                        filter         = 1 ./ (1 + (D_0.^2./D1./D2).^(n_bwf));
                    case char('Gaussian')
                        % Gaussian bandreject filter
                        filter         = 1 - exp(-1/2 * D1.*D2./D_0.^2);
                end
                fft2D_shifted_filter = fft2D_shifted_filter .* filter;
            end
            cf = figure('Name', ['C_scan_amp_2dfftFiltered', '_', F_type, '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(2, 1, 1);
            imagesc(ax, obj.fx*X/nfftx, obj.fy*Y/nffty, 20*log10(abs(fft2D_shifted_filter)/nfftx/nffty)); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %2d fft and filtering
            fft2D_ishift         = ifftshift(fft2D_shifted_filter);
            C_scan_inam_2DFFT    = abs(ifft2(fft2D_ishift));
            % plot 2d fft filtered
            ax     = subplot(2, 1, 2);
            imagesc(Yd, Xd, abs(C_scan_inam_2DFFT(1: m, 1:n)));
            axis image; colormap(jet);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            % display snr
            disp([F_type 'snr and cnr image_filter:']);
            [snr, cnr] = fx_image_snr_2(C_scan_inam_2DFFT(1: m, 1:n), roi);
            disp([snr cnr]);
        end
        
        function apply2dfft_filter(obj, z, PropertyName, filtertype, slicetype)
            % apply 2d fft filter to the C scan by z index.
            % z: z index. unit: data points
            % filtertype: choose the type of the filter:
            % filtertype 1: 'LP', 'HP', 'BP', 'BR', 'norch', 'norchBP'.
            % filtertype 2: 'ideal', 'Butterworth', or 'Gaussian'
            % slicetype: 'z depth', 'surface parallel', 'ply wise'
            % ***
            img = abs(obj.(PropertyName));
            switch slicetype
                case char('z depth')
                    C_scan_inam = img(:, :, z);
                case char('surface parallel')
                    ratio = z;
                    [~, C_scan_inam, ~] = obj.define_parallel_inamCscan(ratio, PropertyName);
                case char('ply wise')
                    layer        = z;
                    inph_ex      = obj.mask_plytrack;
                    [m, n, ~]    = size(img);
                    C_scan_index = NaN(m, n);
                    for i = 1: m
                        for j = 1:n
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
                    C_scan_inam = nan(m, n);
                    for i = 1: m
                        for j = 1:n
                            if obj.rear_I(i, j) > C_scan_index(i, j) && img(i, j, C_scan_index(i, j)) % check the rear surface and assign NaN
                                C_scan_inam(i, j) = inam(i, j, C_scan_index(i, j));
                            else
                                C_scan_inam(i, j) = NaN;
                            end
                        end
                    end
                otherwise % not correct
                    fprint("slicetype shuold be 'z depth', 'surface parallel', or 'ply wise'. ");
            end
            % fillna
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            % plot original
            % display snr
            disp('*************************************************')
            disp('snr and cnr image_origin:')
            % now make 2D fft of original image
            nfftx = 2^nextpow2(m);
            nffty = 2^nextpow2(n);
            nfftx = max(nfftx, nffty); % cannot tackle the issue of different nfftx and nffty.
            nffty = max(nfftx, nffty);
            %
            fft2D         = fft2(C_scan_inam, nfftx, nffty);
            fft2D_shifted = fftshift(fft2D);
            % remove DC
            fft2D_shifted(floor(end/2-4):floor(end/2+6), floor(end/2-4):floor(end/2+6)) = 0;
            X             = -nfftx/2+1: nfftx/2;
            Y             = -nffty/2+1:nffty/2;
            cf = figure('Name', ['2Dspectrum', '_', slicetype, '_', num2str(z)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(1, 1, 1);
            imagesc(X, Y, 20*log10(abs(fft2D_shifted)/m/n)); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %             % *** 2d fft filter
            switch filtertype{1}
                case char('LP')
                    prompt = 'What is the radius of the LP filter? ';
                    D_0 = input(prompt);
                    amp_max_filter = fx_lowpass_2dfft(C_scan_inam, D_0, filtertype{2});
                case char('BP')
                    % Butterworth bandpass filter
                    prompt = 'What is the radius of the BP filter? ';
                    D_0 = input(prompt);
                    prompt = 'What is the width of the BP filter? ';
                    W = input(prompt);
                    prompt = 'What is the order of the Butterworth filter? ';
                    n_bwf = input(prompt);
                    amp_max_filter = fx_bandpass_2dfft(C_scan_inam, D_0, W, filtertype{2}, n_bwf);
                case char('norchBP')
%                     %                     prompt = 'What is the BP points of the BP filter? ';
%                     [B,ind]       = maxk(abs(fft2D_shifted(:)), 10);
%                     [row,col]     = ind2sub(size(fft2D_shifted),ind);
%                     imagesc(abs(fft2D_shifted));
%                     hold on;
%                     scatter(row, col, 15, 'ro');
%                     hold on;
%                     scatter(row, col, 15, 'ro');
                    X_0 = [-31 22 0 22];
                    Y_0 = [0 -22 31 22];
                    prompt = 'What is the order of the Butterworth filter? ';
                    n_bwf = input(prompt);
                    prompt = 'What is the radius of the norchBP filter? ';
                    D_0 = input(prompt);
                    amp_max_filter = fx_norchBP_2dfft(C_scan_inam, X_0, Y_0, D_0, filtertype{2}, n_bwf);
                case char('HP')
                    
                otherwise % no filter
                    amp_max_filter = amp_max;
            end
            % plot 2d fft filtered
            ax     = subplot(2, 1, 1);
            imagesc(ax, Yd, Xd, C_scan_inam);
            axis image;
            colormap(ax, jet);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            ax     = subplot(2, 1, 2);
            imagesc(Yd, Xd, abs(amp_max_filter(1: m, 1:n)));
            axis image; colormap(jet);
            hold on;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
        end
        
        % ******************** 2d fft **********************
        
        function check2dfft_inclphase(obj, z, PropertyName)
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
            set(cf, 'Position', [0, 0, 1200, 400], 'color', 'white');
            ax     = subplot(1, 3, 1);
            Xd     = (1: m) / obj.fx * 1e3;
            Yd     = (1: n) / obj.fy * 1e3;
            imagesc(ax, Yd, Xd, abs(C_scan_inam));
            axis image;
            hold on;
            h = colorbar; colormap(ax, jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y displacement (mm)', 'fontsize', 16);
            set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(gca, 'linewidth', 2);
            ax     = subplot(1, 3, 2);
            log_abs_2dfft                   =  20*log10(abs(fft2D_shifted)/nfftx/nffty);
            log_abs_2dfft(end/2+1, end/2+1) = NaN; % remove DC;
            imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, log_abs_2dfft); % divide by the x and y length: nfftx and nffty
            %             imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, real(fft2D_shifted)/nfftx/nffty); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar; colormap(jet);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            ax     = subplot(1, 3, 3);
            imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, angle(fft2D_shifted)); % divide by the x and y length: nfftx and nffty
            % imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, imag(fft2D_shifted)/nfftx/nffty); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar; colormap(ax, gray);
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Phase (rad.)');
            xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
        end
        
        % ******************** 2D image-wise analysis *******************
        
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
        
        % ********************* 3D ID ********************
        function obj = extract_local_orientation_3DID(obj, PropertyName, wl, orientation, z_theta, z_range, K)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ratios: the ratio determing the index by the front and back surfaces , unit: arb.
            % K: the smoothing applied to the Gabor magnitude responses.
            % sigma: threshould for curvelet denoising, 0 means no application of denoising
            % ***
            % start
            img_temp  = obj.(PropertyName);
            img_temp  = abs(img_temp(:, :, z_range(1): z_range(2)));
            [m, n, l] = size(img_temp); 
%             m_np2     = nextpow2(m);
%             n_np2     = nextpow2(n);    
%             l_np2     = nextpow2(l);
%             % Gaussian filter
%             sigma       = [0.1 0.1 10];
%             img_temp    = imgaussfilt3(img_temp, sigma);
            Inplane_direction           = zeros([m, n, l]);
            Inplane_direction_index_ori = zeros([m, n, l]);
            Inplane_direction_index_zth = zeros([m, n, l]);
            statistic_mean = zeros([m, n, l]);
            statistic_std  = zeros([m, n, l]);
            SFB = 1.5;
            SAR = 0.5;
            %
            disp('convolution calculation...');
            for i = 1:length(orientation)
                for j = 1:length(z_theta)
                    h       = fx_handwritten_3DGabor(orientation(i), wl, z_theta(j), SFB, SAR);
                    %                     fftSize = [2^m_np2, 2^n_np2, 2^l_np2];
                    size_h = size(h);
                    fftSize = size_h + [m n l] - 1;
                    A_k = ifftn( fftn(img_temp, fftSize) .* fftn(real(h), fftSize), 'symmetric' );
                    A_k = imgaussfilt3(A_k, K* 0.5*wl);
                    A_k = A_k(end/2-m/2+1:end/2+m/2, ...
                        end/2-n/2+1:end/2+n/2, ...
                        end/2-l/2+1:end/2+l/2);
                    statistic_mean = statistic_mean + A_k;
                    statistic_std  = statistic_std + A_k.^2;
                    %                     A_k = A_k(wl+1:end-wl, wl+1:end-wl, :);
                    %               A_k = A_k(end/2-m/2+wavelength+1:end/2+m/2-wavelength, ...
                    %                         end/2-n/2+wavelength+1:end/2+n/2-wavelength, ...
                    %                         end/2-l/2+wavelength+1:end/2+l/2-wavelength);
                    % update the maximum
                    compare = A_k>=Inplane_direction; % logic larger or not
                    Inplane_direction_index_ori = Inplane_direction_index_ori.*(1-compare) + i.*compare;
                    Inplane_direction_index_zth = Inplane_direction_index_zth.*(1-compare) + j.*compare;
                    Inplane_direction = Inplane_direction .*(1-compare) + A_k.*compare;
%                     compare = A_k>=Inplane_direction(wl+1:end-wl, wl+1:end-wl, :); % logic larger or not
%                     Inplane_direction_index_ori(wl+1:end-wl, wl+1:end-wl, :)...
%                         = Inplane_direction_index_ori(wl+1:end-wl, wl+1:end-wl, :).*(1-compare) + i.*compare;
%                     Inplane_direction_index_zth(wl+1:end-wl, wl+1:end-wl, :)...
%                         = Inplane_direction_index_zth(wl+1:end-wl, wl+1:end-wl, :).*(1-compare) + j.*compare;
%                     Inplane_direction(wl+1:end-wl, wl+1:end-wl, :)...
%                         = Inplane_direction(wl+1:end-wl, wl+1:end-wl, :).*(1-compare) + A_k.*compare;
                end
                clc;
                fprintf('ID searching progress: %0.2f%%\n',100*i/length(orientation));
            end
            % statistic
            N = length(orientation)*length(z_theta);
            statistic_mean = statistic_mean ./ N;
            statistic_std  = sqrt( (statistic_std - statistic_mean.^2*N)/(N-1) );      
            % remove extra parts
            Inplane_direction_index_ori(1:wl, :, :) = nan;
            Inplane_direction_index_ori(:, 1:wl, :) = nan;
            Inplane_direction_index_ori(end-wl:end, :, :) = nan;
            Inplane_direction_index_ori(:, end-wl:end, :) = nan;
            Inplane_direction_index_zth(1:wl, :, :) = nan;
            Inplane_direction_index_zth(:, 1:wl, :) = nan;
            Inplane_direction_index_zth(end-wl:end, :, :) = nan;
            Inplane_direction_index_zth(:, end-wl:end, :) = nan;
            % statistic removement
            statistic_mean_remove = Inplane_direction<statistic_mean + statistic_std .* 1.5;
            Inplane_direction_index_ori(statistic_mean_remove) = nan;
            Inplane_direction_index_zth(statistic_mean_remove) = nan;
%             [m,n,p] = size(Inplane_direction_index_ori);
%             rear_I_temp = obj.rear_I;
%             front_I_temp = obj.front_I;
%             for i = 1:m
%                 for j =1:n
%                     Inplane_direction_index_ori(i,j, 1:front_I_temp(i,j)) = nan;
%                     if p>rear_I_temp(i,j)
%                         Inplane_direction_index_ori(i,j, rear_I_temp(i,j):p) = nan;
%                     end
%                 end
%             end     
            obj.Inplane_direction_3D_ID = Inplane_direction_index_ori / length(orientation) * 180;
            obj.Inplane_direction_3D_ID_ztheta = Inplane_direction_index_zth / length(orientation) * 180; 
        end
        
        function obj = extract_local_orientation_3DID_plane(obj, PropertyName, wavelength, y_ori, z_ori, z_range, K)
            % extract the local orientation image in each ply
            % PropertyName: 'img', 'img_hil' or 'img_hil_filter' ...
            % wavelength: wavelengths of the filter bank for logGabor filter gabor
            % orientation: orientations of the filter bank for logGabor filter gabor
            % ratios: the ratio determing the index by the front and back surfaces , unit: arb.
            % K: the smoothing applied to the Gabor magnitude responses.
            % sigma: threshould for curvelet denoising, 0 means no application of denoising
            % ***
            % start
            img_temp = obj.(PropertyName);
            img_temp = abs(img_temp(:, :, z_range(1): z_range(2)));
            %             % Gaussian filter
            %             sigma       = [0.1 0.1 10];
            %             img_temp    = imgaussfilt3(img_temp, sigma);
            Inplane_direction = zeros(size(img_temp));
            Inplane_direction_index_ori = zeros(size(img_temp));
            Inplane_direction_index_zth = zeros(size(img_temp));
            SFB = 1.5;
            SAR = 0.5;
            %
            x_theta = 0; % only need to rotate y and z
            disp('convolution calculation...');
            for i = 1:length(y_ori)
                for j = 1:length(z_ori)
                    y_theta = 90 + y_ori(i);
                    z_theta = z_ori(j);
                    thetas = [x_theta, y_theta, z_theta];
                    h = fx_handwritten_3DGabor_plane(thetas, wavelength, SFB, SAR);
                    fftSize = size(img_temp);
                    A_k   = ifftn( fftn(img_temp, fftSize) .* fftn(real(h), fftSize), 'symmetric' );
                    A_k = imgaussfilt3(A_k, K* 0.5*wavelength);
                    % update the maximum
                    compare = A_k>=Inplane_direction; % logic larger or not
                    Inplane_direction_index_ori = Inplane_direction_index_ori.*(1-compare) + i.*compare;
                    Inplane_direction_index_zth = Inplane_direction_index_zth.*(1-compare) + j.*compare;
                    Inplane_direction = Inplane_direction.*(1-compare) + A_k.*compare;
                end
                clc;
                fprintf('ID searching progress: %0.2f%%\n',100*i/length(y_ori));
            end
            %             % remove extra parts
            %             [m,n,p] = size(Inplane_direction_index_ori);
            %             rear_I_temp = obj.rear_I;
            %             front_I_temp = obj.front_I;
            %             for i = 1:m
            %                 for j =1:n
            %                     Inplane_direction_index_ori(i,j, 1:front_I_temp(i,j)) = nan;
            %                     if p>rear_I_temp(i,j)
            %                         Inplane_direction_index_ori(i,j, rear_I_temp(i,j):p) = nan;
            %                     end
            %                 end
            %             end
            obj.Inplane_direction_3D_ID = Inplane_direction_index_ori / length(orientation) * 180;
        end
        
        % ******************** fingerprint of ply *********************
        function obj = fingerprint_2dfft_fitsurface(obj, ratios, PropertyName)
            % find 2dfft features along z depth,
            % fitting the surface profiles
            % start
            [lxc, lyc, lzc] = size(obj.(PropertyName));
            fft2d_abs       = NaN(lxc, lyc, lzc);
            fft2d_pha       = NaN(lxc, lyc, lzc);
            fft2d_mask      = NaN(lxc, lyc, lzc, 2);
            index_container = NaN(lxc, lyc);
            % incl. location, phase, amplitude ...
            phase_info = nan(length(ratios), 5);
            row = [70;82;63;89;94;58;63;89;76;76];
            col = [89;63;63;89;76;76;81;71;57;95];
            ind = sub2ind(size(fft2d_pha), row, col); % convert to linear index
            for i = 1:length(ratios)
                ratio = ratios(i);
                %                 ratio = ratios(max(1, i-3): min(i+3, end)); % the ratio to determine the index
                if i~=length(ratios)
                    ratio_next = ratios(i+1);
                end
                index_container_prev = index_container;
                [~, C_scan_inam_para, index_container] = obj.define_parallel_inamCscan(ratio, PropertyName);
                [~, ~, index_container_next]           = obj.define_parallel_inamCscan(ratio_next, PropertyName);
                % no filter
                C_scan_inam_para_denoise = C_scan_inam_para;
                % 2d fft
                % now make 2D fft of original image
                %                 nfftx         = 2^nextpow2(lxc);
                %                 nffty         = 2^nextpow2(lyc);
                %
                fft2D         = fft2(C_scan_inam_para_denoise, lxc, lyc);
                fft2D_shifted = fftshift(fft2D);
                % normalize
                fft2D_shifted(floor(end/2+1), floor(end/2+1)) = nan;
                fft2D_shifted = fft2D_shifted/max(abs(fft2D_shifted), [], 'all', 'omitnan');
                %
                if i==1
                    upper_index_bound = index_container;
                    lower_index_bound = round(index_container/2 ...
                        +index_container_next/2);
                elseif i==length(ratios)
                    upper_index_bound = round(index_container/2 ...
                        + index_container_prev/2);
                    lower_index_bound = index_container;
                else
                    upper_index_bound = round(index_container/2 ...
                        + index_container_prev/2);
                    lower_index_bound = round(index_container/2 ...
                        +index_container_next/2);
                end
                % mask
                fft2d_abs_oneplane = 20*log10(abs(fft2D_shifted));
                fft2d_pha_oneplane = angle(fft2D_shifted);
                %
                mask = fft2d_abs_oneplane<mean(fft2d_abs_oneplane(:), 'omitnan')+6;
                fft2d_mask_abs_oneplane = fft2d_abs_oneplane;
                fft2d_mask_abs_oneplane(mask) = nan;
                fft2d_mask_pha_oneplane = fft2d_pha_oneplane;
                fft2d_mask_pha_oneplane(mask) = nan;
                for ii = 1:lxc
                    for jj = 1:lyc
                        fft2d_abs(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj))...
                            = fft2d_abs_oneplane(ii, jj);
                        fft2d_pha(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj))...
                            = fft2d_pha_oneplane(ii, jj);
                        fft2d_mask(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj))...
                            = fft2d_mask_abs_oneplane(ii, jj); % abs
                        fft2d_mask(ii, jj, upper_index_bound(ii, jj): lower_index_bound(ii, jj))...
                            = fft2d_mask_pha_oneplane(ii, jj);     % pha
                    end
                end
                %
                cos_pha = cos(fft2d_pha_oneplane);
                cos_pha = imgaussfilt(cos_pha);
                phase_info(i, :) = cos_pha(ind);
                % end loop in one Cscan
                clc;
                disp(['Procedure: ', num2str(round(100*ratios(i)/ratios(end))), '%']);
            end
            %
            disp(['number of nan ', num2str(round(sum(isnan(fft2d_abs), 'all')))]);
            
            % average abs
            ave_fft2d_abs = mean(fft2d_abs, 3, 'omitnan');
            [B,ind]       = maxk(ave_fft2d_abs(:), 10);
            [row,col]     = ind2sub(size(ave_fft2d_abs),ind);
            imagesc(ave_fft2d_abs);
            hold on;
            scatter(row, col, 15, 'ro');
            hold on;
            scatter(row, col, 15, 'ro');
            
            %
            obj.fft2d_abs_fingerprint = fft2d_abs;
            obj.fft2d_pha_fingerprint = fft2d_pha;
            obj.fft2d_mask_fingerprint = fft2d_mask;
        end
        
        function obj = fingerprint_2dfft_fitsurface_v2(obj, ratios, PropertyName)
            % find 2dfft features along z depth,
            % fitting the surface profiles
            % start
            [lxc, lyc, lzc] = size(obj.(PropertyName));
            fft2d_abs       = NaN(lxc, lyc, length(ratios));
            fft2d_pha       = NaN(lxc, lyc, length(ratios));
            fft2d_mask      = NaN(lxc, lyc, length(ratios), 2);
            % incl. location, phase, amplitude ...
            row = [70;82;63;89;94;58;63;89;76;76];
            col = [89;63;63;89;76;76;81;71;57;95];
            phase_info = nan(length(ratios), length(row));
            ind = sub2ind(size(fft2d_pha), row, col); % convert to linear index
            for i = 1:length(ratios)
                ratio = ratios(i);
                [~, C_scan_inam_para, ~] = obj.define_parallel_inamCscan(ratio, PropertyName);
                % no filter
                C_scan_inam_para_denoise = C_scan_inam_para;
                % ****** fit a plane
                [x_mesh, y_mesh]         = meshgrid(1:lxc, 1:lyc);
                plane_surfacefit         = fit([x_mesh(:), y_mesh(:)], C_scan_inam_para_denoise(:), 'poly11');
                plane_surface            = plane_surfacefit(x_mesh, y_mesh);
                C_scan_inam_para_denoise = C_scan_inam_para_denoise - plane_surface.';
                % 2d fft
                % now make 2D fft of original image
                %                 nfftx         = 2^nextpow2(lxc);
                %                 nffty         = 2^nextpow2(lyc);
                %
                fft2D         = fft2(C_scan_inam_para_denoise, lxc, lyc);
                fft2D_shifted = fftshift(fft2D);
                % normalize
                fft2D_shifted = fft2D_shifted/lxc/lyc;
%                 fft2D_shifted(floor(end/2+1), floor(end/2+1)) = nan;
%                 fft2D_shifted = fft2D_shifted/max(abs(fft2D_shifted), [], 'all', 'omitnan');
                % mask
                fft2d_abs_oneplane = 20*log10(abs(fft2D_shifted));
                fft2d_pha_oneplane = angle(fft2D_shifted);
                %
                mask = fft2d_abs_oneplane<mean(fft2d_abs_oneplane(:), 'omitnan')+6;
                fft2d_mask_abs_oneplane = fft2d_abs_oneplane;
                fft2d_mask_abs_oneplane(mask) = nan;
                fft2d_mask_pha_oneplane = fft2d_pha_oneplane;
                fft2d_mask_pha_oneplane(mask) = nan;
                for ii = 1:lxc
                    for jj = 1:lyc
                        fft2d_abs(ii, jj, i) = fft2d_abs_oneplane(ii, jj);
                        fft2d_pha(ii, jj, i) = fft2d_pha_oneplane(ii, jj);
                        fft2d_mask(ii, jj, i, 1) = fft2d_mask_abs_oneplane(ii, jj); % abs
                        fft2d_mask(ii, jj, i, 2) = fft2d_mask_pha_oneplane(ii, jj);     % pha
                    end
                end
                %
                cos_pha = cos(fft2d_pha_oneplane);
                cos_pha = imgaussfilt(cos_pha);
                phase_info(i, :) = cos_pha(ind);
                % end loop in one Cscan
                clc;
                disp(['Procedure: ', num2str(round(100*ratios(i)/ratios(end))), '%']);
            end
            %
            disp(['number of nan ', num2str(round(sum(isnan(fft2d_abs), 'all')))]);
  
            % average abs
            ave_fft2d_abs = mean(fft2d_abs, 3, 'omitnan');
            [~,ind]       = maxk(ave_fft2d_abs(:), 10);
            [row,col]     = ind2sub(size(ave_fft2d_abs),ind);
            imagesc(ave_fft2d_abs);
            hold on;
            scatter(row, col, 45, 'ro');
            %
            obj.fft2d_abs_fingerprint = fft2d_abs;
            obj.fft2d_pha_fingerprint = fft2d_pha;
            obj.fft2d_mask_fingerprint = fft2d_mask;
        end
        
        function [best_object, best_location] = fingerprint_oneply(obj, z, PropertyName)
            % check fingerpint of one woven ply
            % z: z index. unit: data points
            % PropertyName: property to sue
            % ***
            % fillna
            %             C_scan_inam = abs(obj.(PropertyName)(:, :, z));
            [~, C_scan_inam, ~] = obj.define_parallel_inamCscan(z, PropertyName);
            C_scan_inam(isnan(C_scan_inam)) = mean(C_scan_inam, 'all', 'omitnan');
            [m, n] = size(C_scan_inam);
            %
            % now make 2D fft of original image
            %             nfftx         = 2^nextpow2(m);
            %             nffty         = 2^nextpow2(n);
            % cannot solve the problem when nfftx!=nffty !
            %             nfftx         = max(nfftx, nffty);
            %             nffty         = max(nfftx, nffty);
            %
            %             fft2D         = fft2(C_scan_inam, nfftx, nffty);
            %             fft2D_shifted = fftshift(fft2D);
            %             % Gaussian Filter Response Calculation
            %             X              = -nfftx/2+1: nfftx/2;
            %             Y              = -nffty/2+1:nffty/2;
            % ****** fit a plane
            [x_mesh, y_mesh]    = meshgrid(1:m, 1:n);
            plane_surfacefit    = fit([x_mesh(:), y_mesh(:)], C_scan_inam(:), 'poly11');
            plane_surface       = plane_surfacefit(x_mesh, y_mesh);
            C_scan_inam_detrend = C_scan_inam - plane_surface.';
            %
            h = max(C_scan_inam_detrend(:))/2 + abs(min(C_scan_inam_detrend(:))/2);
            omega = 2 * pi * 600 / obj.fx;
            best_object   = 1e10;
            best_location = [1 1];
            % ******
            pha1_range = 30;
            pha2_range = 30;
            %             % *** debug
%             surface = fx_SinuSurface(h, omega, best_location(1)*2*pi/pha1_range, ...
%                 best_location(2)*2*pi/pha2_range, x_mesh, y_mesh);
%             figure, surf(C_scan_inam_detrend);
%             colormap(jet);
%             hold on;
%             surf(surface,'FaceAlpha',0.5, 'EdgeColor','none');%             object_fit = nan(pha1_range, pha2_range);
            for i = 1:pha1_range
                for j = 1:pha2_range
                    phase1 = 2*pi * i/pha1_range;
                    phase2 = 2*pi * j/pha2_range;
                    surface = fx_SinuSurface(h, omega, phase1, phase2, x_mesh, y_mesh);
                    object_t = sum((C_scan_inam_detrend - surface).^2, 'all', 'omitnan');
                    % % save
                    %                     object_fit(i,j) = object_t;
                    % % search for best solution
                    if best_object > object_t
                        best_location = [i j];
                        best_object = object_t;
                    end
                end
            end
            % ***********
            %             % **** plot original
            %             cf = figure('Name', ['C_scan_amp' , '_', num2str(z)]);
            %             set(cf, 'Position', [0, 0, 1200, 400], 'color', 'white');
            %             ax     = subplot(1, 3, 1);
            %             Xd     = (1: m) / obj.fx * 1e3;
            %             Yd     = (1: n) / obj.fy * 1e3;
            %             imagesc(ax, Yd, Xd, abs(C_scan_inam));
            %             axis image;
            %             hold on;
            %             h = colorbar; colormap(ax, jet);
            %             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            %             xlabel('\fontname {times new roman} X displacement (mm) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Y displacement (mm) ', 'fontsize', 16);
            %             set(gca, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(gca, 'linewidth', 2);
            %             ax     = subplot(1, 3, 2);
            %             log_abs_2dfft                   =  20*log10(abs(fft2D_shifted)/nfftx/nffty);
            %             log_abs_2dfft(end/2+1, end/2+1) = NaN; % remove DC;
            %             imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, log_abs_2dfft); % divide by the x and y length: nfftx and nffty
            %             %             imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, real(fft2D_shifted)/nfftx/nffty); % divide by the x and y length: nfftx and nffty
            %             axis image;
            %             h = colorbar; colormap(jet);
            %             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            %             xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            %             set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(ax, 'linewidth', 2);
            %             ax     = subplot(1, 3, 3);
            %             imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, angle(fft2D_shifted)); % divide by the x and y length: nfftx and nffty
            %             % imagesc(obj.fx*X/nfftx, obj.fy*Y/nffty, imag(fft2D_shifted)/nfftx/nffty); % divide by the x and y length: nfftx and nffty
            %             axis image;
            %             h = colorbar; colormap(ax, gray);
            %             set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Phase (rad.)');
            %             xlabel('\fontname {times new roman} X wavenumber (1/m) ', 'fontsize', 16);
            %             ylabel('\fontname {times new roman} Y wavenumber (1/m)', 'fontsize', 16);
            %             set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            %             set(ax, 'linewidth', 2);
        end
        
        % ************************ 3D fft *********************************
        function obj = nDfft(obj, ra, rb, rc, F_type)
            % perform 3d fft on the original img
            img_temp = obj.img_hil; % apply fft on the analytic-signal directly
            %             [lx, ly, lz] = size(img_temp);
            %             l = nextpow2(lx);
            %             m = nextpow2(ly);
            %             n = nextpow2(lz);
            %             ndfft = fftn(img_temp, [2^l 2^m 2^n]);
            %             ndfft = fftshift(ndfft);
            img_filtered2 = fx_lowpass_3dfft(img_temp, ra, rb, rc, F_type);
            
            % f = obj.fs * (0:(l / 2)) / l;
            obj.img_3dfft  = ndfft;
            obj.img_fftfilter =img_filtered2;
        end
        
        function obj = hilbert_property(obj, propertyName)
            img_temp    = obj.(propertyName);
            img_hilbert = nan(size(img_temp));
            for i = 1: size(img_hilbert, 1)
                %                 for j = 1:size(img_hilbert, 2)
                %
                %                 end
                bscan     = squeeze(img_temp(i, :, :));
                bscan_hil = hilbert(bscan.'); %　hilbert for each colume at 2D matrix
                img_hilbert(i, :, :) = bscan_hil.'; %　hilbert for each colume at 2D matrix
                clc;
                disp([num2str(i) '/' num2str(size(img_hilbert, 1))]);
            end
            obj.img_hil_filter = img_hilbert;
        end
        
        % ************************ determine defect size ******************
        function obj = amplitude_drop_method(obj, PropertyName, zrange, drop, filtertype)
            % determine thoverall size of image damagas
            % PropertyName: the name of the property in the object
            % zrange: the z axis range to find the max amplitudes
            % drop: the ratio of drop, -3 dB, -6 dB or .....
            % filtertype: type of the filter.
            img_temp = abs(obj.(PropertyName));
            img_temp = img_temp(:,:,zrange(1):zrange(2)); % select the range
            amp_max  = max(img_temp, [], 3, 'omitnan');
            % plot the 2d spectrum
            [m, n] = size(amp_max);
            
            %             % ************* debug and EDA ***************** %
            %             inam_ex  = abs(obj.img_hil);
            %             inph_ex  = angle(obj.img_hil);
            %             infq_ex  = diff(unwrap(inph_ex, [], 3), 1, 3) /2 / pi * obj.fs;
            %
            %             infq_ex_single = infq_ex>10e6;
            %             infq_ex_single(infq_ex>10e6) = 100;
            %             infq_ex_single(infq_ex<2e6) = -100;
            %             infq_ex_single(infq_ex>2e6 & infq_ex<10e6) = 0;
            % %             volumeViewer(inph_ex);
            %             close all;
            %             temp   = inph_ex(:,:,600); % select the range
            %             Xd     = (1: m) / obj.fx * 1e3;
            %             Yd     = (1: n) / obj.fy * 1e3;
            %             figure, imagesc(Yd, Xd, temp);
            %             % now make 2D fft of original image
            %             nfftx = 2^nextpow2(m);
            %             nffty = 2^nextpow2(n);
            %             nfftx = max(nfftx, nffty); % cannot tackle the issue of different nfftx and nffty.
            %             nffty = max(nfftx, nffty);
            %             %
            %             fft2D         = fft2(temp, nfftx, nffty);
            %             fft2D_shifted = fftshift(fft2D);
            %             X              = -nfftx/2+1: nfftx/2;
            %             Y              = -nffty/2+1:nffty/2;
            %             figure, imagesc(X, Y, 20*log10(abs(fft2D_shifted)/m/n)); % divide by the x and y length: nfftx and nffty
            %             axis image;
            %             % ************* end EDA ***************** %
            
            % now make 2D fft of original image
            nfftx = 2^nextpow2(m);
            nffty = 2^nextpow2(n);
            nfftx = max(nfftx, nffty); % cannot tackle the issue of different nfftx and nffty.
            nffty = max(nfftx, nffty);
            %
            fft2D         = fft2(amp_max, nfftx, nffty);
            fft2D_shifted = fftshift(fft2D);
            X              = -nfftx/2+1: nfftx/2;
            Y              = -nffty/2+1:nffty/2;
            cf = figure('Name', ['2Dspectrum', '_', num2str(zrange)]);
            set(cf, 'Position', [0, 0, 600, 800], 'color', 'white');
            ax     = subplot(1, 1, 1);
            imagesc(X, Y, 20*log10(abs(fft2D_shifted)/m/n)); % divide by the x and y length: nfftx and nffty
            axis image;
            h = colorbar;
            set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (dB)');
            set(ax, 'Fontname', 'times new Roman', 'Fontsize', 16);
            set(ax, 'linewidth', 2);
            %             % *** 2d fft filter
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
            
            %             [Y, X]       = size(amp_max);
            %             cw           = [5 10 20 40];
            %             filtStruct   = createMonogenicFilters(Y, X, cw, 'lg', 0.66);
            %             [m1, m2, m3] = monogenicSignal(amp_max, filtStruct);
            % %             [m1, m2, m3] = monogenicSignal_nofilter(amp_max, filtStruct);
            %             close all;
            %             % Local energy (calculated on a per-scale basis)
            %             LE = localEnergy(m1, m2, m3);
            %             % Local phase (calculated on a per-scale basis)
            %             LP = localPhase(m1, m2, m3);
            %             % Local orientation (calculated on a per-scale basis)
            %             % Only need to pass the odd parts (m2,m3) as even part (m1) is irrelevant
            %             LO = localOrientation(m2, m3);
            %
            %             for i = 1:size(LE, 4)
            %                 figure();
            %                 C_scan_inam_denoise = sqrt(LE(:, :, :, i)); % local amplitdue = sqrt(local energy)!
            %                 imagesc(C_scan_inam_denoise), axis image;
            %                 colormap jet;
            %                 h = colorbar;
            %                 set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');
            %             end
            % find the average peak amplitudes
            amp_max_mean = mean(amp_max_filter, 'all');
            log_amp_max  = 20*log10(amp_max_filter);
            groudtruth   = log_amp_max <= 20*log10(amp_max_mean) + drop;
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
        
        
    end
end

