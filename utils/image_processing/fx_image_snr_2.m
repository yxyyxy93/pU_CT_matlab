function [snr, cnr] = fx_image_snr_2(image, roi)
% definition of cnr from DOI: 10.1109/TUFFC.2019.2956855
% definition of cnr from DOI: 10.1109/TUFFC.2019.2956855
% roi: the region of interest, 1st row: x, 2nd row: y.

xv = roi(1,1):roi(1,2);
yv = roi(2,1):roi(2,2);

ROI = image(xv, yv); % 
ROO = image;
ROO(xv, yv) = nan;


mu_i = mean(ROI.^2, 'all', 'omitnan');
mu_o = mean(ROO.^2, 'all', 'omitnan');

snr = 20*log10(mean(ROI, 'all', 'omitnan') / std(ROO, 0, 'all', 'omitnan'));

% cnr is defined as CNR = ΔS/σ, where ΔS is the signal enhancement and σ
% is the noise standard deviation [11].
% note that here A_max is used to approxmate the ΔS!

cnr = 10*log10(mu_i / mu_o); 

end

