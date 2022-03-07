function img_filtered2 = fx_lowpass_3dfft(img, ra, rb, rc, F_type)
% lowpass filter using 2d fft
% img: the image to process
% ra, rb, rc: radius of the filter
% F_type: type of the filter: 'ideal', 'Butterworth', or 'Gaussian'
% n_bwf: order of the Butterworth filter

% now make 3D fft of original image
[lx, ly, lz] = size(img);
nfftx = 2^nextpow2(lx);
nffty = 2^nextpow2(ly);
nfftz = 2^nextpow2(lz);
ndfft = fftn(img, [nfftx nffty nfftz]);
ndfft = fftshift(ndfft);
            
% cannot solve the problem when nfftx!=nffty !
nfftx         = max(nfftx, nffty);
nffty         = max(nfftx, nffty);

%
X = -nfftx/2+1: nfftx/2;
Y = -nffty/2+1:nffty/2;
Z = -nfftz/2+1:nfftz/2;

[Xmesh, Ymesh, Zmesh] = meshgrid(X, Y, Z);
% design the filter
D = sqrt(Xmesh.^2/ra^2 + Ymesh.^2/rb^2 + Zmesh.^2/rc^2);

switch F_type
    case char('ideal')
        Lo = exp(-D.^2);
%         Lo = exp(-D.^2)./(2 * R_Lo).^2);
    case char('sinc')
        % Butterworth bandreject filter
        ;     
end
%
%3d fft and filtering
fft3D_shifted_filter = ndfft .* Lo;

% % volumeViewer(abs(fft3D_shifted_filter), );
% s = orthosliceViewer(abs(fft3D_shifted_filter), 'Colormap', jet, 'DisplayRange', [0 .5]...
%     , 'ScaleFactors', [1 1 0.2]);
% % for debug
% figure, imagesc(20*log10(abs(fft2D_shifted_filter))), axis image;

%             fft2D_shifted_filter = fft2D_shifted;
fft3D_ishift         = ifftshift(fft3D_shifted_filter);
img_filtered1        = ifftn(fft3D_ishift, 'symmetric');
img_filtered2        = img_filtered1(1:lx, 1:ly, 1:lz);

end

