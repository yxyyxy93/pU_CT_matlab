function [B, bw, f_lower, f_upper]= fx_bandwidth_octaves(sigma0, f0)
% calculate the bandwidth from the sigma_0 governing the bandwidth in log Gabor filter.

B       = sqrt(8 / log(2)) * abs(log(sigma0));
ratio   = 2^B;

% by formula
f_lower = exp(log(f0) - 0.5*B*log(2));
f_upper = f_lower * ratio;

bw      = (f_upper - f_lower) / f0;

disp([num2str(sigma0) ':']);
disp(f_lower);
disp(f_upper);
disp(bw);

% by direct calculation
B       = exp((2*log(2))^0.5 * log(sigma0));
f_lower = f0 * B;
f_upper = f0 * B^(-1);
bw      = (f_upper - f_lower) / f0;

end

