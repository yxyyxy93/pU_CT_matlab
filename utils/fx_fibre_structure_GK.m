% This program generates the fibre structure 2D gaussian filter.
% there is some problem with the rotation function,
% so angle should be 0 right now.

function kernel = fx_fibre_structure_GK(n, s, angle)

x = -1/2:1/(n-1):1/2;
[~, X] = meshgrid(x,x);
kernel = exp( -(X.^2+X.^2)/(2*s^2) );

kernel = imrotate(kernel, angle);

kernel = kernel - mean(kernel(:));

end
