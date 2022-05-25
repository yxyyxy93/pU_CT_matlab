function y = gausswin2XY(Nx,Ny, alpha)
% function y = gausswin2XY(Nx,Ny, alpha)
% 
%  function generates 2D Gaussian window using the following mathematical
%  expression 
% y = exp(-1/2*(R/stdev).^2);
%
% where 
% R = sqrt(NX.^2+NY.^2); 
% stdev = (Nx-1)/(2*alpha);
% 
% Nx and Ny is the number of points for the window in x and y direction 
% 
%           (c) ambrozin@agh.edu.pl
% 
% example: 
% Nx = 25;      Ny = 25; 
% alph_f = 1; 
% y = gausswin2XY(Nx,Ny, alph_f); 
% figure(1),surf(y)
% title(['\alpha = ' num2str(alph_f) ])


nx = -(Nx-1)/2:(Nx-1)/2;
ny = -(Ny-1)/2:(Ny-1)/2;
% nx = nx./max(nx); 
ny = ny./max(ny).*max(nx);
[NX NY] = meshgrid(nx,ny); 
% NY = meshgrid(ny); 

R = sqrt(NX.^2+NY.^2); 

stdev = (Nx-1)/(2*alpha);
y = exp(-1/2*(R/stdev).^2);

