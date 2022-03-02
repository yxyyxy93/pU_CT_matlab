function [snr, cnr] = fx_image_snr_1(image)
% from wikipedia 
% definition of SNR is as the reciprocal of the coefficient of variation, 
% i.e., the ratio of mean to standard deviation of a signal or measurement:
% notice that such an alternative definition is only useful for variables 
% that are always non-negative (such as photon counts and luminance)

snr = mean(image(:)) / std(image(:));

% cnr is defined as CNR = ΔS/σ, where ΔS is the signal enhancement and σ
% is the noise standard deviation [11].
% note that here A_max is used to approxmate the ΔS!

cnr = max(image(:))/ std(image(:)); 

end

