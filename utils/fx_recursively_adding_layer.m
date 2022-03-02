function [R_fwd, T_back, R_back, T_fwd] = fx_recursively_adding_layer(R_12, R_21, R_23, R_32, T_12, T_21, T_23, T_32, alpha, k, l)
%**** xiaoyu yang
%**** checked on 30/10/2019

% resin-rich inter-ply layers are initially converted into single 
% interfaces with complex reflection and transmission coefficients and then
% each ply, plus its following inter-ply interface, is successively 
% combined into overall laminate reflection (R) and 
% transmission (T) coefficients (for acoustic pressure).

% k = 2*pi/lambda; % wavenumber of the layer

% * there may be some mistakes in the paper, '-' should be '+' in the denominator!!
R_fwd  = R_12 + R_23.*T_12.*T_21.*exp(-2.*l.*(alpha+1i.*k)) ./ (1-R_21.*R_23.*exp(-2.*l.*(alpha+1i.*k)));
R_back = R_32 + R_21.*T_32.*T_23.*exp(-2.*l.*(alpha+1i.*k)) ./ (1-R_23.*R_21.*exp(-2.*l.*(alpha+1i.*k)));

% * there may be some mistakes in the paper, there should be no exp in the molecule!!
T_fwd  = T_12.*T_23.*exp(-l.*(alpha+1i.*k)) ./ (1-R_23.*R_21.*exp(-2.*l.*(alpha+1i.*k)));
T_back = T_32.*T_21.*exp(-l.*(alpha+1i.*k)) ./ (1-R_23.*R_21.*exp(-2.*l.*(alpha+1i.*k)));

%********from PhD thesis appendix
% R_fwd  = R_12 - R_12.*T_12.*T_21.*exp(-2.*l.*(alpha+1i.*k)) ./ (1-R_12.*R_12.*exp(-2.*l.*(alpha+1i.*k)));

end

