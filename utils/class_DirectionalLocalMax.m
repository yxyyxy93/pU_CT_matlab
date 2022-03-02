classdef class_DirectionalLocalMax < handle
    %UNTITLED3 Summary of this class goes here
    % linear_window (pixel_number, angle)
    % angle from (-pi/2 pi/2)
  
    properties
        pixelsnum
        angle
        in
        out
        processed
    end
    
    methods

        function DLM = class_DirectionalLocalMax(in, linear_window)
            DLM.pixelsnum = linear_window(1);
            DLM.angle = linear_window(2);
            DLM.in = in;
            DLM.out = uint8(zeros(size(DLM.in)));
            DLM.processed = false(size(DLM.in));
        end

        % change the linear window
        function change_LWin(DLM, linear_window)
            DLM.pixelsnum = linear_window(1);
            DLM.angle = linear_window(2);
        end
        
        function process_singlepixel(DLM, ii, jj)
            %UNTITLED3 Construct an instance of this class
            %   Detailed explanation goes here
            
            if DLM.processed(ii, jj)==0 %false
                
                %for each pixel in linear window
                for win_index = (1-DLM.pixelsnum)/2:(DLM.pixelsnum-1)/2
                    % not compare itself
                    
                    if (win_index==0)
                        continue;
                    end
                    
                    if (abs(DLM.angle)<=pi/4) % symmetry
                        x_win = ii + win_index;
                        y_win = jj + round(win_index*tan(DLM.angle));
                    else
                        x_win = ii + round(win_index/tan(DLM.angle));
                        y_win = jj + win_index;
                    end
                    
                    if (DLM.in(ii,jj)<=DLM.in(x_win, y_win))
                        return;
                    else
                        DLM.processed(x_win, y_win) = 1; % set true
                    end
                    
                end  
                
%                 DLM.out(ii, jj) = DLM.in(ii, jj);
                DLM.out(ii, jj) = 1;
            end
            
        end
        
        function process_image(DLM)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
       
            % for each pixel in image 
            for ii = ceil(DLM.pixelsnum/2):size(DLM.in,1)-ceil(DLM.pixelsnum/2)
                for jj = ceil(DLM.pixelsnum/2):size(DLM.in,2)-ceil(DLM.pixelsnum/2)
                    DLM.process_singlepixel(ii, jj);
                end
            end
            
        end
        
        function show_window(DLM, ii, jj, fignum)
            linear_window_x = zeros(1, DLM.pixelsnum);
            linear_window_y = zeros(1, DLM.pixelsnum);
            
            for i=1:DLM.pixelsnum
                win_index = (1-DLM.pixelsnum)/2+i-1;
                if (abs(DLM.angle)<=pi/4) % symmetry
                    linear_window_x(i) = ii + win_index;
                    linear_window_y(i) = jj + round(win_index*tan(DLM.angle));
                else
                    linear_window_x(i) = ii + round(win_index/tan(DLM.angle));
                    linear_window_y(i) = jj + win_index;
                end
            end
            figure(fignum);
            scatter(linear_window_x, linear_window_y);
        end
        
        % reset the out
        function reset_out(DLM)
            DLM.out = uint8(zeros(size(DLM.in)));
            DLM.processed = false(size(DLM.in));
        end
    end
        
end

