function [R_12, R_21, T_12, T_21] = fx_Pressure_R_T(s1, s2, theta1)
%********** xiaoyu yang *************
%************** 30/10/2019 **********

%*** s1, Z1 represent the fisrt layer where the incidence is on.
switch nargin
    case 2
        %*********** acoustic impedance ***********
        Z_1 = s1.rho1*s1.cl;
        Z_2 = s2.rho1*s2.cl;
    case 3
        Z_1 = s1.rho1*s1.cl / cosd(theta1);
        theta2 = asin(sin(theta1) / s1.cl * s2.cl);
        Z_2 = s2.rho1*s2.cl / cosd(theta2);
end
        
%*********** acoustic impedance ***********
Z_1 = s1.rho1*s1.cl;
Z_2 = s2.rho1*s2.cl;

%*******Pressure Reflection and Transmission Coeff*******
R_12 = (Z_2-Z_1) / (Z_2+Z_1);
R_21 = (Z_1-Z_2) / (Z_2+Z_1);
T_12 = 2*Z_2 / (Z_2+Z_1);
T_21 = 2*Z_1 / (Z_2+Z_1);

end

