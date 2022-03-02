function [R, theta_l, theta_t]=fx_R_theta_omega(c0,cl,ct,alpha_l,alpha_t,warray,h1,rho0,rho1,theta)
%theta 
%20170918
%20171102 
%20181107 
%20181107 

theta_l=asin(sin(theta).*cl./c0);
theta_t=asin(sin(theta).*ct./c0);

% constant
costhetal=cos(theta_l);
costhetat=cos(theta_t);
sinthetal=sin(theta_l);
sinthetat=sin(theta_t);
tanthetat=tan(theta_t);

Lw=length(warray);
R=zeros(1,Lw);
% warray

kl=warray./cl.*(1+1i.*alpha_l);
kt=warray./ct.*(1+1i.*alpha_t);
%º
klz=kl.*costhetal;
ktz=kt.*costhetat;

%a21=i.*(2.*sinthetat.^2.*sin(klz.*h1)./tanthetal-tanthetat.*cos(2.*theta_t).*sin(ktz.*h1));
a21=1i.*(2.*sinthetat.*costhetal.*sin(klz.*h1).*ct./cl-tanthetat.*cos(2.*theta_t).*sin(ktz.*h1));
a22=cos(2.*theta_t).*cos(klz.*h1)+2.*sinthetat.^2.*cos(ktz.*h1);
a23=(costhetal.*sin(klz.*h1)+tanthetat.*sinthetal.*sin(ktz.*h1))./(warray.*rho1.*cl);
a31=-2.*1i.*warray.*rho1.*ct.*sinthetat.*cos(2.*theta_t).*(cos(ktz.*h1)-cos(klz.*h1));
a32=-warray.*rho1.*(cl.*cos(2.*theta_t).^2./costhetal.*sin(klz.*h1)+4.*ct.*sinthetat.^2.*costhetat.*sin(ktz.*h1));
a41=-warray.*rho1.*ct.^2.*(4.*costhetal.*sinthetat.^2.*sin(klz.*h1)./cl+cos(2.*theta_t).^2.*sin(ktz.*h1)./costhetat./ct);
  


R=(a32.*a41-a31.^2+(warray.*rho0.*c0./cos(theta)).^2.*(a23.*a41-a21.^2))./(a32.*a41-a31.^2-(warray.*rho0.*c0./cos(theta)).^2.*(a23.*a41-a21.^2)-2.*1i.*(warray.*rho0.*c0./cos(theta)).*a41.*a22-a21.*a31);

end

