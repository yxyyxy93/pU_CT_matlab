function [R, anguler_1D, max_angle_I, xp, theta] = fx_Slice_RD_1DA(inam_C_scan_mask)
% radon transform and reduce the 2D RDT iamge to 1D angular distribution
% select the max point in 1D angular distribution as well

% inam_C_scan_mask is already the prepared slice by "fx_extractSlice" 

theta = 0:180;
[R, xp] = radon(inam_C_scan_mask, theta);

% calculating the gradient by central-different kernel
[Gx, ~] = imgradientxy(R, 'intermediate');

anguler_1D = sum(abs(Gx), 1);
[~, max_angle_I] = max(anguler_1D);

end

