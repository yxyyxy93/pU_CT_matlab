function [wl, ratio]= fx_kernelsize_R2wl(R, SAR, SFB)
% input the radius of a disp-shaped area 
% output the wavelength of the 2D Gabor filter with the same kernel size
% 
% sigma_x = sqrt(R.^2*SAR);

ratio = sqrt(SAR/2) / (1/pi * sqrt(log(2)/2) * (2^SFB + 1)/(2^SFB - 1));

wl =  R * ratio;

end

