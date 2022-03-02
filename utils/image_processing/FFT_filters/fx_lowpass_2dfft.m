function img_filtered2 = fx_lowpass_2dfft(img, R_Lo, F_type)
% lowpass filter using 2d fft
% img: the image to process
% R_Lo: radius of the filter
% F_type: type of the filter: 'ideal', 'Butterworth', or 'Gaussian'
% n_bwf: order of the Butterworth filter

[m, n] = size(img);
% now make 2D fft of original image
nfftx         = 2^nextpow2(m);
nffty         = 2^nextpow2(n);
% cannot solve the problem when nfftx!=nffty !
nfftx         = max(nfftx, nffty);
nffty         = max(nfftx, nffty);
%
fft2D         = fft2(img, nfftx, nffty);
fft2D_shifted = fftshift(fft2D);
X              = -nfftx/2+1: nfftx/2;
Y              = -nffty/2+1:nffty/2;
[Xmesh, Ymesh] = meshgrid(Y, X);
% design the filter
D              = sqrt(Xmesh.^2 + Ymesh.^2);
switch F_type
    case char('ideal')
        Lo             = exp(-(D.^2)./(2 * R_Lo).^2);
    case char('sinc')
        % Butterworth bandreject filter
        ;     
end
%
%2d fft and filtering
fft2D_shifted_filter = fft2D_shifted .* Lo;

% for debug
% figure, imagesc(20*log10(abs(fft2D_shifted_filter))), axis image;

%             fft2D_shifted_filter = fft2D_shifted;
fft2D_ishift         = ifftshift(fft2D_shifted_filter);
img_filtered1        = ifft2(fft2D_ishift, 'symmetric');
img_filtered2        = img_filtered1(1:m, 1:n);

end

