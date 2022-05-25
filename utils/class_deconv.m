classdef class_deconv
    % used for deconvolution
    % including the Autoregressive
    
    properties
        recorded_signal
        ref_signal
        fs
    end
    
    methods
        function obj = class_deconv(recorded_signal, ref_signal, fs)
            % a class to implement the deconvolution, including the
            % AutoRregressive extrapolation
            obj.recorded_signal = recorded_signal;
            obj.ref_signal = ref_signal;
            obj.fs = fs;
        end
        
        function deconv = deconvolution(obj)
            % normal deconvolution
            % just divide the spectrum, ignoring the noises
            n = 2^nextpow2(length(obj.recorded_signal));
            deconv_spectrum = fft(obj.recorded_signal, n) ./ (fft(obj.ref_signal, n));
            deconv = ifft(deconv_spectrum);
        end
        
        function [deconv, deconv_spectrum] = wiener_deconlution(obj, Q_factor, fft_padding)
            % wiener deconvolution
            % Q_factor: the Q factor in wiener deconvolution function
            % output the deconv_spectrum for AR extrapolation
            if nargin == 1
                Q_factor = 1e-2; % should be 10^(-2) here
            end
            kernel = obj.ref_signal;
            ori_signal = obj.recorded_signal;
            shift = round(length(kernel) / 2);
            n = 2^nextpow2(length(obj.recorded_signal)) * 2^fft_padding;
            H = fft(kernel, n);
            Y_omega  = fft(ori_signal, n);
            % This constant is sometimes called the ‘‘noise desensitizing factor’’
            Q = sqrt(Q_factor * max(H .* conj(H)).^2);
            deconv_spectrum = Y_omega .* conj(H) ./ (H .* conj(H) + Q.^2);
            devolved = ifft(deconv_spectrum);
            % The shift maybe not needed
            deconv = circshift(devolved, shift);
        end
        
        function deconv = wienerdeconv_ARextrap(obj, Q_factor, f_windows, k, method, fft_padding)
            % Autoregressive extrapolation
            % Q_factor: Q_factor in wiener deconvolution
            % f_windows: the spectrum windows to fit the model; [number of the windows, 2(f1, f2)]
            % k: order the the AR model
            % method: 'bg' or 'yk'
            % calculate the wiener deconvolved specturm
            if nargin == 1
                Q_factor = 1e-2; % should be 10^(-2) here
            end  
            kernel = obj.ref_signal;
%             % move the max to the center
%             [~, maxI] = max(kernel);
%             kernel = kernel(1: 2* maxI);
            ori_signal = obj.recorded_signal;
            shift      = round(length(kernel) / 2);
            n          = fft_padding;
            H          = fft(kernel, n);
            Y_omega    = fft(real(ori_signal), n);
            % This constant is sometimes called the ‘‘noise desensitizing factor’’
            Q                           = sqrt(Q_factor * max(H .* conj(H)).^2);
            deconv_spectrum  = Y_omega .* conj(H) ./ (H .* conj(H) + Q.^2);
            %  
            ns                           = length(deconv_spectrum);
            AR_spetrum_ave   = 0;          
            for i = 1:size(f_windows, 1)
                f_window              =  f_windows(i, :);
                f_window_index   = round(ns * f_window / obj.fs); % transfer the frequency window to the index.
                m                           = f_window_index(1);
                n                            = f_window_index(2);
                % for debug
                %             disp(f_window_index);
                spectrum_window = deconv_spectrum(m: n);
                % debug
                %             figure, plot(abs(deconv_spectrum));
                %             hold on;
                %             scatter(f_window_index, abs(deconv_spectrum(f_window_index)));
                %             title("The deconvolved spectrum");
                %             % AR model
                % choose methods
                if strcmp(method, 'bg')
                    [a, ~, ~] = arburg(spectrum_window, k);
                elseif strcmp(method, 'yk')
                    [a, ~, ~] = aryule(spectrum_window,k);
                else
                    error("please choose method as 'bg' or 'yk'. ");
                end
                % transverse a and rc
                a      = a';
                a_conj = conj(a);
                % extrapolate the miss values outside the window
                AR_spetrum       = zeros(ns,1);
                AR_spetrum(m: n) = deconv_spectrum(m: n);
                % p = 1, 2, ... m - 1
                for p = m-1: -1: 1
                    AR_spetrum(p) = - sum( AR_spetrum(p + 1: p + k) .* a(2: end));
                end
                % q = n + 1, n + 2, ... Ns
                for q = n + 1: 1: ns
                    AR_spetrum(q) = - sum( AR_spetrum(q - 1: -1: q - k) .* a_conj(2: end) ) ;
                end
                %             % Reconstruct full spectrum, above Nyquist frequency
                %             Yx = zeros(2*length(AR_spetrum)-1,1);
                %             Yx(1:length(AR_spetrum))= AR_spetrum(1:end);
                %             Yx(length(AR_spetrum)+1:end)= conj(flip(Yx(2:length(AR_spetrum))));
                %             % debug
                %             figure, plot(abs(AR_spetrum)) ;
                %             %
                AR_spetrum_ave = AR_spetrum_ave + AR_spetrum;
            end
            AR_spetrum_ave = AR_spetrum_ave / size(f_windows, 1);
            devolved       = real(ifft(AR_spetrum_ave));
            deconv         = circshift(devolved, shift);
        end
    end
end

