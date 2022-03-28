function img_filtered2 = fx_norchBP_2dfft(img, Xv_0, Yv_0, D_0, F_type, n_bwf)
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
X              = -nfftx/2+1: nfftx/2;
Y              = -nffty/2+1:nffty/2;
[Xmesh, Ymesh] = meshgrid(Y, X);
%
fft2D         = fft2(img, nfftx, nffty);
fft2D_shifted = fftshift(fft2D);
fft2D_shifted_filter = fft2D_shifted;

for i = 1:length(Xv_0)
    X_0 = Xv_0(i);
    Y_0 = Yv_0(i);
    D1  = sqrt((Xmesh-X_0).^2 + (Ymesh-Y_0).^2);
    D2  = sqrt((Xmesh+X_0).^2 + (Ymesh+Y_0).^2); % symmetry
    switch F_type
        case char('ideal')
            filter = ones(nfftx, nffty);
            filter(D1>D_0 & D2>D_0) = 0;
        case char('Butterworth')
            % Butterworth bandreject filter
            filter = 1 - 1 ./ (1 + (D_0.^2./D1./D2).^(n_bwf));
        case char('Gaussian')
            % Gaussian bandreject filter
            filter = exp(-1/2 * D1.*D2./D_0.^2);
    end
    %
    %2d fft and filtering
    fft2D_shifted_filter = fft2D_shifted_filter .* filter;
end
% for debug
% figure, imagesc(20*log10(abs(fft2D_shifted_filter))), axis image;

%             fft2D_shifted_filter = fft2D_shifted;
fft2D_ishift  = ifftshift(fft2D_shifted_filter);
img_filtered1 = ifft2(fft2D_ishift, 'symmetric');
img_filtered2 = img_filtered1(1:m, 1:n);

end

