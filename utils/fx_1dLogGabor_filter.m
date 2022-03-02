function G = fx_1dLogGabor_filter(f, f0, sigma)
% ADDME generate a Log-Gabor filter 
% f: the freq. space
% f0 and sigma  are the parameters of the filter 
% f0 gives the center frequency of the filter.
% sigma  affects the bandwidth of the filter.

mol = ( log(abs(f) ./ f0) ).^2;
den = 2 * ( log(sigma) ).^2;

G = exp(- mol ./ den);

end

