function [alpha, A_intial] = fx_calcualate_attcoef(A, A0, h, z1, z2)
% calculate attenuation coefficient by A=A0*exp(-alpha*z)

T12 = 2*z2/(z2+z1);
R12 = (z2-z1)/(z2+z1);
R21 = (z1-z2)/(z2+z1);
T21 = 2*z1/(z2+z1);

k     = A/A0;
alpha = log(abs(R12*k/(T12*R21*T21)))/(-h);

A_intial = A0 * abs((T12*R21*T21) / R12);

end

