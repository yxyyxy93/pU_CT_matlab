function im_noise_remov = ACARPENOS(im , size_patches, logNFAthresh, reg, roi)
% code for the paper 
% "An a-contrario approach to quasi-periodic noise removal"
% Frederic Sur
% IEEE International Conference on Image Processing (ICIP)
% 2015
%
% usage:
% im: input image (double)
% size_patches: size of the patches (L in the paper)
% logNFAthresh: detection threshold on log10(NFA) (set to 0 in the paper)
% reg (optional argument): 0= no TV regularization, 1= constrained TV regularization

disp(' ')
disp('AC-ARPENOS v. 1.2')
disp('Please cite:')
disp('Frederic Sur')
disp('An a-contrario approach to quasi-periodic noise removal')
disp('IEEE International Conference on Image Processing (ICIP)')
disp('2015.')
disp(' ')
disp('Wait a few seconds...')
drawnow('update')

switch nargin
    case 3
        reg=0;
end
        
% main constants (cf. paper)
step_patches=round(size_patches/2); % translation between patches (in pixels)
f2=8/size_patches; 


% display image power spectrum (log scale) 
Fim=fft2(im);
Specim=abs(fftshift(Fim)).^2;
figure; imagesc(log(Specim)); colormap(gray), colorbar
title('power spectrum of the initial image (log scale)')
colorbarspectrum=[min(min(2*log(abs(fftshift(Fim))))), max(max(2*log(abs(fftshift(Fim)))))];
[rows, cols] = size(im);


% calculate minimum power spectrum
P=numel(1:step_patches:(rows-size_patches))*numel(1:step_patches:(cols-size_patches));  % total number of patches
disp(' ')
H=hanning(size_patches)*hanning(size_patches)';  % patches are weighted by a Hann window to remove border effects
logspeclist=zeros(size_patches,size_patches,P);
n=0;
for i=1:step_patches:(rows-size_patches)
  for j=1:step_patches:(cols-size_patches)
  n=n+1;
  imagette=H.*im(i:i+size_patches-1,j:j+size_patches-1);
  logspeclist(:,:,n)=2*log(abs(fftshift(fft2(imagette))));
  end
end
logspecmin=min(logspeclist,[],3);

% display minimum power spectrum 
figure; imagesc(logspecmin); colormap(gray), colorbar
title('minimum power spectrum on patches (log scale)')


% calculate empirical background law and NFA
kstep=floor(size_patches/8);
[fx, fy] = meshgrid([-size_patches/2:size_patches/2-1]/size_patches, [-size_patches/2:size_patches/2-1]/size_patches);
Freq=sqrt(fx.^2+fy.^2); 
f1=max(Freq(:));
proba=zeros(size(imagette));
NumElem=zeros(size(imagette));
for k=0:kstep-1
        J=find((Freq>=f2+k*(f1-f2)/kstep)&(Freq<=f2+(k+1)*(f1-f2)/kstep)); % define the rings of the paper
        numelJ=numel(J);
        K=find((Freq>=f2+0.75*k*(f1-f2)/kstep)&(Freq<=f2+1.25*(k+1)*(f1-f2)/kstep));  % Fourier coefficients are cumulated from a slightly larger region than J
        numelK=numel(K); 
        NumElem(J)=numelJ;
        Background=zeros(1,numelK*P);
        for p=1:P
            logspecimagette=logspeclist(:,:,p);
            Background(1+(p-1)*numelK:p*numelK)=[reshape(logspecimagette(K),[1,numelK])];
        end
        for j=1:numel(J)
            proba(J(j))=sum(Background>=logspecmin(J(j)))/(numelK*P);
        end
end

logNFA=log10(kstep)+log10(NumElem)+P*log10(proba);

% display NFA map
map=[ 1:-.01:.51, zeros(1,10);  1:-.01:.51, zeros(1,10); zeros(1,50),1:-.05:0.55]';
figure, imagesc(logNFA,[-50,10]), colorbar, colormap(map)
title('logarithm of the NFA (base 10)')

logNFA(Freq<f2)=10;
disp(['Most meaningful coefficient:  log_10(NFA) = ',num2str(min(logNFA(:)))])
disp(' ')


% Outlier map 
Peaks=exp(logspecmin).*(logNFA<logNFAthresh);
Peaks(Freq<f2)=0;  % NFA is not defined for low-frequency components

% Interpolation of the outlier map to the original image dimension
Peaks_b=zeros(size(Peaks)+2);
Peaks_b(2:end-1,2:end-1)=Peaks;
[fx, fy] = meshgrid([-size_patches/2-1:size_patches/2]/size_patches, [-size_patches/2-1:size_patches/2]/size_patches);
[x, y] = meshgrid([-cols/2:cols/2-1]/cols, [-rows/2:rows/2-1]/rows);
Out=interp2(fx,fy,Peaks_b,x,y);

% Notch filter design
H = fspecial('gaussian',[12 12],2);
Notch_filter=sym_conv2(double(Out>0),H);
Notch_filter= fftshift(fft2(real(ifft2(ifftshift(Notch_filter))))); % (symmetrization to get rid of rounding errors)
M=max(Notch_filter(:));
if (M>0) 
    Notch_filter=Notch_filter/M; 
end


% Method noise estimation
Noise=real(ifft2(Fim.*ifftshift(Notch_filter)));


% Final display
figure; imagesc(2*log(abs(fftshift(Fim).*(1-Notch_filter))),colorbarspectrum); colormap(gray), colorbar;
title('corrected power spectrum')
figure; imagesc(im); colormap(jet), colorbar
title('original image')
im_noise_remov=im-Noise;
figure; imagesc(im_noise_remov); colormap(jet), colorbar;
title('denoised image')
figure; imagesc(Noise), colormap(jet); colorbar
title('noise component')

%
disp('*************************************************')
disp('snr and cnr image_ACARPENOS:')
[snr, cnr] = fx_image_snr_2(im_noise_remov, roi);
disp([snr cnr]);

% Constrained TV regularization
if (reg)
    disp('Constrained TV regularization')
    disp(' ')
    niter=500; % max. iterations
    epsilon=0.1; % regularization parameter;
    Ima=im_noise_remov; % steepest descent initialization
    Gstep=1;
    Gx=diff(Ima([end 1:end],:),1,1); 
    Gy=diff(Ima(:,[end 1:end]),1,2); 
    regnorm=sqrt(epsilon^2+Gx.^2+Gy.^2);
    TV=mean(regnorm(:));
    disp(['Mean Total Variation before regularization: TV = ',num2str(TV)])
    for i=1:niter
        tempx=Gx./regnorm; tempy=Gy./regnorm;
        Gxx=diff(tempx([1:end 1],:),1,1); 
        Gyy=diff(tempy(:,[1:end 1]),1,2); 
        divergence=Gxx+Gyy;
        GradProj=-ifft2(fft2(divergence).*ifftshift(Notch_filter),'symmetric'); % projected gradient
        Ima2=Ima-Gstep*GradProj; % constant step size
        Gx2=diff(Ima2([end 1:end],:),1,1); 
        Gy2=diff(Ima2(:,[end 1:end]),1,2); 
        regnorm2=sqrt(epsilon^2+Gx2.^2+Gy2.^2);
        TV2=mean(regnorm2(:)); % new value of the mean TV
        while (TV2>TV) 
            Gstep=Gstep/2; 
            disp([' Gstep = ',num2str(Gstep)]); 
            Ima2=Ima-Gstep*GradProj; 
            Gx2=diff(Ima2([end 1:end],:),1,1); 
            Gy2=diff(Ima2(:,[end 1:end]),1,2); 
            regnorm2=sqrt(epsilon^2+Gx2.^2+Gy2.^2);
            TV2=mean(regnorm2(:)); 
        end
        if (TV-TV2<1e-4)
            break;
        else
            TV=TV2;
            regnorm=regnorm2;
            Ima=Ima2;
            Gx=Gx2;
            Gy=Gy2;
        end
    end
    disp(['Mean Total Variation after ',num2str(i),' iterations of gradient descent: TV = ',num2str(TV)])
    disp(' ')
    im_noise_remov_regul=Ima;
    figure; imagesc(2*log(abs(fftshift(fft2(Ima)))),colorbarspectrum); colormap(gray), colorbar;
    title('power spectrum of the regularized denoised image')
    figure; imagesc(im_noise_remov_regul,[0 255]); colormap(gray), colorbar;
    title('regularized denoised image')
    figure; imagesc(im-im_noise_remov_regul), colormap(gray); colorbar
    title('Difference between original image and regularized denoised image')
    figure; imagesc(im_noise_remov-im_noise_remov_regul), colormap(gray); colorbar
    title('Difference between denoised images before and after regularization')
    
    im_noise_remov=im_noise_remov_regul;

end
