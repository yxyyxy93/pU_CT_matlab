function img_filtered2 = fx_bandpass_2dfft(img, D_0, W, F_type, n_bwf)
% bandpass filter using 2d fft
% img: the image to process
% D_0: radius of the filter
% W: width of the filter 
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
        filter = ones(nfftx, nffty);
        filter(D>(D_0+W/2) | D<(D_0-W/2)) = 0;      
    case char('Butterworth')
        % Butterworth bandreject filter
        filter         = 1 - 1 ./ (1 + (D*W./(D.^2-D_0^2)).^(2*n_bwf));
    case char('Gaussian')
        % Gaussian bandreject filter
        filter         = exp(-1/2 * ((D.^2-D_0^2)./D./W).^2);
end
%
%2d fft and filtering
fft2D_shifted_filter = fft2D_shifted .* filter;

% for debug
% figure, imagesc(20*log10(abs(fft2D_shifted_filter))), axis image;

%             fft2D_shifted_filter = fft2D_shifted;
fft2D_ishift         = ifftshift(fft2D_shifted_filter);
img_filtered1        = ifft2(fft2D_ishift, 'symmetric');
img_filtered2        = img_filtered1(1:m, 1:n);

end

