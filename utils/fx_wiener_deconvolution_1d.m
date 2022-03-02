function devolved = fx_wiener_deconvolution_1d(ori_signal, kernel, q_factor)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
if nargin == 2
    q_factor = 1e-2;
end

shift = round(length(kernel) / 2);
kernel = [kernel, zeros(1, length(ori_signal) - length(kernel))];  % zero pad the kernel to same length
H = fft(kernel) / 1; % shrink the kernel
Y_omega  = fft(ori_signal);
% This constant is sometimes called the ‘‘noise desensitizing factor’’
Q = sqrt(q_factor * max(H .* conj(H)).^2);
devolved = real(ifft(Y_omega .* conj(H) ./ (H .* conj(H) + Q.^2)));
% devolved = real(ifft(fft(ori_signal) ./ (H);
devolved = circshift(devolved, shift);

end