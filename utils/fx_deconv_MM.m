function [x, cost] = fx_deconv_MM(y, phi, wfun, H, Nit)
% [x, cost] = deconv_MM(y, phi, wfun, H, Nit)
% Sparse deconvolution
% Cost function : 0.5 * sum(abs(y-H(x)).^2) + sum(phi(x));
%
% INPUT
% y - noisy data
% H - banded convolution matrix [sparse data structure]
% phi - cost function
% wfun - function x / phi’(x) - needed for quadratic MM
% Nit - number of iterations
%
% OUTPUT
% x - deconvolved sparse signal
% cost - cost function history
% Reference: Penalty and Shrinkage Functions for Sparse Signal Processing
% Ivan Selesnick, selesi@poly.edu, 2012
% http://eeweb.poly.edu/iselesni/lecture_notes/
% Algorithm: majorization-minimization with banded system solver.
y          =  y(:); % convert to column vector
cost       =  zeros(1, Nit); % cost function history
[M, N]     = size(H);
x          = y; % initialization
g          = H' * y; % H’*y
for k = 1:Nit
    Lam      = spdiags(wfun(x), 0, N, N); % Lam : diagonal matrix
    F          = speye(M) + H*Lam*H'; % F : banded matrix
    x          = Lam * (g - (H'*(F\(H*(Lam*g))))); % update x (solve banded system)
    cost(k) = 0.5*sum(abs(y-H*x).^2) + sum(phi(x)); % cost function value
end

end