function m = fx_sinc2d(xymax, nrc, freq)
%Create the 2-D array 
%Input: [xymax,nrc,freq]
% Output: [m] an array of size(nrc,nrc) representing the sin(x)/x function

m(row,col) = sin(pi*R*freq)/(pi*R*freq);

% where R is the distance from the center of the array. 
% The array is [nrc,nrc] with X(1,1)=Y(1,1)= - xymax and X(nrc,nrc)=Y(nrc,nrc)=xymax.


end

