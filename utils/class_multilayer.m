classdef class_multilayer
    % a multilayer structure in an Ascan
    %   Detailed explanation goes here
    
    properties
        structure
        n_layers
        Rc
        spikes
    end
    %    properties (Dependent)
    %       Modulus
    %    end
    methods
        function specimen = class_multilayer(structure)
            specimen.structure = structure;
            specimen.n_layers = length(structure);
            specimen.Rc = 0;
        end
        
        % calculate reflectionCoefficiency
        function specimen = reflectionCoefficiency_oblique(specimen, upper_media, back_media, freq, theta)
            % define the count of the layer
            global count_layer;
            count_layer = 1;
            
            % calculate A
            c0   = upper_media.cl;
            rho0 = upper_media.rho1;
            [A, ~] = fx_R_theta_omega_all(c0, freq * 2 * pi, rho0, theta, specimen.structure{1});
 
            % remember
            M = containers.Map();
            
            for i = 1: specimen.n_layers - 1        
%                 h = waitbar(0, ['layer',  num2str(count_layer), ': Please wait...']);
                if (isKey(M, specimen.structure{i + 1}.Material))
                    a2 = M(specimen.structure{i + 1}.Material);
                else
                    [a2,  ~] = fx_R_theta_omega_all(c0, freq * 2 * pi, rho0, theta, specimen.structure{i + 1});
                    M(specimen.structure{i + 1}.Material) = a2;
                    disp(specimen.structure{i + 1}.Material);
                end
                
                for j = 1:length (freq)
                    A(j, :, :) = squeeze(A(j, :, :)) * squeeze(a2(j, :, :));
%                     if mod(j, 1000)==0
%                         disp([num2str(j), '/', num2str(length(freq))]);
%                     end
                end
%                 waitbar((i - 1) / specimen.n_layers - 1, h);    
%                 close(h);
%                 count_layer = count_layer + 1;
            end
            
            a32 = squeeze(A(:, 3, 2));
            a41 = squeeze(A(:, 4, 1));
            a23 = squeeze(A(:, 2, 3));
            a21 = squeeze(A(:, 2, 1));
            a31 = squeeze(A(:, 3, 1));
            a22 = squeeze(A(:, 2, 2));
            w = freq * 2 * pi;
            w = transpose(w);
         
            R_fwd = (a32.*a41 - a31.^2 + (w.*rho0.*c0./cos(theta)).^2.*(a23.*a41- a21.^2)) ...
                ./(a32.*a41 - a31.^2 - (w.*rho0.*c0./cos(theta)).^2.*(a23.*a41 - a21.^2) - ...
                2.*1i.*(w.*rho0.*c0./cos(theta)).*a41.*a22 - a21.*a31);
            
            % conjugate and transpose if needed
            specimen.Rc = R_fwd';
        end
        
        % calculate reflectionCoefficiency
        function specimen = reflectionCoefficiency(specimen, upper_media, back_media, freq)
            %************ one resin layer, fisrt layer***********
            s_fisrt                  = specimen.structure{1};
            [R_12, R_21, T_12, T_21] = fx_Pressure_R_T(upper_media, s_fisrt);  % fisrt boundary from infinium space
            %********* adding layers recusively ****************
            for i=1:specimen.n_layers-1
                s_add                          = specimen.structure{i};
                s_back                         = specimen.structure{i+1};
                [R_23, R_32, T_23, T_32]       = fx_Pressure_R_T(s_add, s_back);
                l                              = s_add.h1;
                k                              = 2.*pi.*freq./s_add.cl;
                db_att                         = s_add.alpha_l;
                [R_fwd, T_back, R_back, T_fwd] = fx_recursively_adding_layer(R_12, R_21, R_23,R_32, T_12, T_21, T_23, T_32, db_att.*freq, k, l);
                %*******recurse the R and T **************
                R_12 = R_fwd;
                R_21 = R_back;
                T_12 = T_fwd;
                T_21 = T_back;
            end
            
            %************ one resin layer, last layer***********
            s_last = specimen.structure{specimen.n_layers};
            [R_23, R_32, T_23, T_32] = fx_Pressure_R_T(s_last, back_media);  % s_last to back_media is the last boundary to infinium space
            l = s_last.h1;
            k = 2.*pi.*freq./s_last.cl;
            db_att = s_last.alpha_l;
            [R_fwd, ~, ~, ~] = fx_recursively_adding_layer(R_12, R_21, R_23,R_32, T_12, T_21, T_23, T_32, db_att.*freq, k, l);       
            specimen.Rc = R_fwd;
        end
        
        % wave incidence
        function echo = reflected_wave(specimen, pulse)
            S1_a = ifft(fft(pulse, length(specimen.Rc)).*specimen.Rc); % amplitude and phase change
            echo = S1_a;
        end
        
        function loc = structrue2loc(specimen)
            loc = zeros(1,specimen.n_layers);
            
            layer = specimen.structure{1};
            loc(1) = layer.h1/layer.cl;
            
            for i=2:specimen.n_layers
                layer = specimen.structure{i};
                loc(i) = loc(i-1) + layer.h1/layer.cl;
            end   
        end
        
    end
    
end