classdef class_reference_signal
    % read and process a reference signal
    properties
        filename
        Props
        setup
        oimg_uint8
        fs
    end
    
    methods
        function obj = class_reference_signal(filename)
            % filename: the path and the filename of the  
            my_tdms_struct = TDMS_getStruct(filename);
            obj.filename = filename;
            [obj.Props, obj.setup] = fx_read_setup(my_tdms_struct);
            % load the sampling frequency
            obj.fs = str2double(obj.setup.SampleRateHz{1});
            % read the signal 
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
                        obj.oimg_uint8 = zeros(numel(fieldnames(my_tdms_struct))-1,...
                            numel(fieldnames(index_struct))-1, ...
                            length(scan_point_struct.data));
                        obj.oimg_uint8(i-1, j-2, :) = scan_point_struct.data;
                    else
                        obj.oimg_uint8(i-1, j-2, :) = scan_point_struct.data;
                    end
                end
            end
            
        end
        
        function img_shape = show_img_shape(obj)
            % display the shape to the oimg
            img_shape = size(obj.oimg_uint8);
            disp(img_shape);
        end
        
        function show_hilbert_ref(obj, x, y)
            % display the decomposed analytical signal of A scan
            % x, y: index of the A scan
            % reuse the fx_showAS
            ascan     = squeeze(obj.oimg_uint8(x, y, :));
            ascan_hil = hilbert(ascan);
            fx_showAS((1:length(ascan)) / obj.fs, ascan_hil, obj.fs);
        end
    end
end

