function [infq, inph, inam, hs]= fx_ht_inFqPhAm_filter(s, fs, f0, sigma)
% This code is designed to calculate the instantaneous frequency (IF), phase and amplitude using Hilbert Transform.
% Input:
% s: signal (either the row or the column matrix will be ok)
% fs: sampling frequency (Hz)
% Output:
% infq: IF (Hz), inph: IP (rad), inam: IA (arb.)
%------------------------
% written by Linshan Jia (jialinshan123@126.com)
% Xi'an Jiaotong University (XJTU)
% 2018-10-26
%
% modified by Xiaoyu Yang (xiaoyuyang@zju.edu.cn)
% Ghent University & Zhejiang University
% 2019-11-29
%------------------------

if length(size(s))==1 || size(s, 1)==1 || size(s, 2)==1
    hs = hilbert(s);
elseif (length(size(s))==3)
    % some problem when using hilbert() to 3D matrix directly!!
    s_size = size(s);
    hs     = zeros(s_size);
    for i=1:s_size(1)
        for j=1:s_size(2)
            Ascan = squeeze(s(i, j, :));
            % calculate the AS in frequency domain
            Fo = fft(Ascan);
            % transverse Fo for dot multiply
            Fo = Fo';
            fx_shift = fftshift(Fo);
            n = length(Ascan);
            fshift = (-n/2: n/2-1) * (fs / n); % zero-centered frequency range
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
            hs(i, j,:) = Ascan_as;
        end
    end
else
    msg = 'The input must be an 1D array or a 3D matrix.';
    error(msg);
end

inph  = atan2(imag(hs), real(hs));
phase = unwrap(inph);
dp    = diff(phase);
inam  = abs(hs);

if length(size(s))==1 || size(s, 1)==1 || size(s, 2)==1 % 1D array
    infq = ([0, dp]+[dp, 0])./2;
    infq = infq./2./pi.*fs;
    infq = [infq(2), infq(2:end-1), infq(end-1)];
elseif (length(size(s))==3) % 3D matrix
    infq = (cat (3, zeros( size(dp(:, :, 1))), dp) + cat (3, dp, zeros( size(dp(:, :, 1)))))./2;
    infq = infq./2./pi.*fs;
    infq = cat(3, infq(:, :, 2), infq(:, :, 2:end-1), infq(:, :, end-1));
end
% trim
% inam = inam(1: length(s));
% inph = inph(1: length(s));
% infq = infq(1: length(s));

end