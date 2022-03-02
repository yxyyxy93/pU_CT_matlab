function [a, R] = fx_R_theta_omega_all(c0,w, rho0, theta, s)
%theta 
%20170918-R
%20170926

cl = s.cl;
ct = s.ct;
alpha_l = s.alpha_l;
alpha_t = s.alpha_t;
rho1 = s.rho1;
h1 = s.h1;

theta_l=asin(sin(theta).*cl / c0);
theta_t=asin(sin(theta).*ct / c0);

%constant
costhetal = cos(theta_l);
costhetat = cos(theta_t);
sinthetal = sin(theta_l);
sinthetat = sin(theta_t);
tanthetat = tan(theta_t);
tanthetal = tan(theta_l);
% k 
% the unit of alpha is /m/Hz, here unit of w is rad/s = 2 * pi  / s = 2 * pi Hz
% also the formula should be modified 

kl = w / cl + 1i.* w / 2 / pi * alpha_l;
kt = w / ct + 1i.* w / 2 / pi * alpha_t;

klz = kl.*cos(theta_l);
ktz = kt.*cos(theta_t);

%a21=i.*(2.*sin(theta_t)^2.*sin(klz.*h1)/tan(theta_l)-tan(theta_t).*cos(2.*theta_t).*sin(ktz.*h1));
a21 = 1i.*(2.*sinthetat.*costhetal.*sin(klz.*h1).*ct./cl-tanthetat.*cos(2.*theta_t).*sin(ktz.*h1));
a22 = cos(2.*theta_t).*cos(klz.*h1)+2.*sinthetat.^2.*cos(ktz.*h1);
a23 = (costhetal.*sin(klz.*h1)+tanthetat.*sinthetal.*sin(ktz.*h1))./(w.*rho1.*cl);
a31 = -2.*1i.*w.*rho1.*ct.*sinthetat.*cos(2.*theta_t).*(cos(ktz.*h1)-cos(klz.*h1));
a32 = -w.*rho1.*(cl.*cos(2.*theta_t).^2./costhetal.*sin(klz.*h1)+4.*ct.*sinthetat.^2.*costhetat.*sin(ktz.*h1));
a41 =-w.*rho1.*ct.*(4.*costhetal.*sinthetat.^2.*sin(klz.*h1).*ct./cl+cos(2.*theta_t).^2.*sin(ktz.*h1)./costhetat);
    
a11=2 .* sinthetat.^2.*cos(klz.*h1) + cos(2.*theta_t).*cos(ktz.*h1);
a12=1i .* (tanthetal .* cos(2.*theta_t) .*sin(klz.*h1) - sin(2.*theta_t) .* sin(ktz.*h1));
a13=1i.* sinthetal .* (cos(ktz.*h1) - cos(klz.*h1)) ./ (w .* rho1 .* cl);
a14=(tanthetal .* sinthetat .* sin(klz.*h1) + costhetat .* sin(ktz.*h1)) ./ (w .* rho1 .* ct);
    
a24=a13;
a33=a22;
a34=a12;
a42=a31;
a43=a21;
a44=a11;
    
% a1 = cat(3, a11, a12, a13, a14);
% a2 = cat(3, a21, a22, a23, a24);
% a3 = cat(3, a31, a32, a33, a34);
% a4 = cat(3, a41, a42, a43, a44);
% 
% a = squeeze(cat(4, a1, a2, a3, a4));
a = ones(length(w), 4, 4);
for i = 1:length(w)
    a(i, :, :) = [
        a11(i) a12(i) a13(i) a14(i);
        a21(i) a22(i) a23(i) a24(i);
        a31(i) a32(i) a33(i) a34(i); 
        a41(i) a42(i) a43(i) a44(i)
        ];
end

R = (a32.*a41 - a31.^2 + (w.*rho0.*c0./cos(theta)).^2.*(a23.*a41-a21.^2)) ...
    ./(a32.*a41 - a31.^2 - (w.*rho0.*c0./cos(theta)).^2.*(a23.*a41 - a21.^2) - ...
    2.*1i.*(w.*rho0.*c0./cos(theta)).*a41.*a22 - a21.*a31);

end

