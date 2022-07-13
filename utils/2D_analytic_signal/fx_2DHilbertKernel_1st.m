function kernel = fx_2DHilbertKernel_1st(x, y, s)
% The first order 2D convolution kernels will be calculated by

ss = s .^ 2;
kk = x .^ 2 + y .^ 2;

kernel = 1 ./ (2 .* pi .* (ss + kk) .^ 1.5);

end

