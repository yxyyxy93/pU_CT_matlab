function kernel = fx_2DHilbertKernel_2nd(x, y, s)
% The first order 2D convolution kernels will be calculated by

ss = s .^ 2;
kk = x .^ 2 + y .^ 2;

d = kk .^ 2 .* (ss + kk).^1.5 .* 2 .* pi;

kernel = -(s .* (2 .* ss + 3 .* kk) - 2 .* (ss + kk) .^ 1.5) ./ d;
kernel(d==0) = 0;

end

