classdef class_process_IRTdata < class_process_RFdata
    % sub class of class_process_RFdata
    properties
        c_l % the linear anisotropy metrics
        c_s % the solid anisotropy metrics
    end
    
    methods
        function obj = class_process_IRTdata(filename)
            %
            obj = obj@class_process_RFdata(filename);
        end
        
        function obj = read_origin_data_IRT(obj)
            fns = fieldnames(obj.data);        
            img_origin   = obj.data.(fns{1});
            [lx, ly, lz] = size(img_origin);
            img_hilbert  = zeros(lx, ly, lz);
            for i = 1:lx
                for j = 1:ly
                    img_hilbert(i, j, :) = hilbert(img_origin(i, j, :));
                end
                disp(i);
            end
            obj.img = single(img_origin);
            obj.img_hil = single(img_hilbert);
            % temporary assigned
            obj.fx = 1;
            obj.fy = 1;
            obj.fs = 1;
        end
        
        % structure tensor process
        function obj = structural_tensor_IRT(obj, sigma1, sigma2, PropertyName, ds_rate)
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
            inph = angle(obj.(PropertyName));
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
            [cp, cs, cl] = fx_decomSS_checkfeature(ST_in);
            %             [~, anglex, angley, anglez] = fx_decomSS_surface(ST_in, obj.front_I, obj.rear_I);
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
            obj.c_p     = cell2mat(cp);
            obj.c_s = cell2mat(cs);
            obj.c_l = cell2mat(cl);
            %             obj.angle_z = cell2mat(anglez);
            % time record
            timeElapsed = toc;
            disp(['convert to arrays: ', num2str(timeElapsed)]);
        end
    end
end

