function [Filter] = gauss_win_filter(fmin, fmax,f2tof1, fs, Nsampl,varargin)
% ____filter parameters__________
% fmin=100e3;
% fmax5=5e6;              fmax10=10e6;
% deltaf=1/2048/5e-9;     Nmin=fmin/deltaf;
% Nmax5=fmax5/deltaf;     Nmax10=fmax10/deltaf;
% f2tof1=1.1;
% ____end filter parameters__________
%
%
%
%
% Ex (1) 
% fs = 200e6; fmin = 0; fmax = 10e6; Nsampl = 2048; f2tof1 = 1; 
% fv = (0:Nsampl-1) .* fs/Nsampl; 
% [F] = gass_win_filter(fmin, fmax,f2tof1, fs, Nsampl);
% figure(3), plot(fv, F)
%
% Ex (2) 
% fs = 200e6; fmin = 1e6; fmax = 10e6; Nsampl = 2048; f2tof1 = 1.2; 
% fv = (0:Nsampl-1) .* fs/Nsampl; 
% [F] = gass_win_filter(fmin, fmax,f2tof1, fs, Nsampl);
% figure(3), plot(fv, F)

if isempty(varargin)
   nPow = 2; 
else 
    nPow = varargin{1};
end
% nPow=4;
deltaf=1/Nsampl.*fs;        Nmin=fmin/deltaf*4;           Nmax10=fmax/deltaf;

Filter = zeros(1, Nsampl);
if fmin == 0
    for i=0:(Nsampl/2),    Filter(i+1)=(exp(-(i/Nmax10)^2));  end
    for i=(Nsampl/2) : (Nsampl-1), Filter(i+1)=(exp(-((Nsampl-i)/Nmax10)^2));  end
else  
    for i=0:(Nsampl/2),    Filter(i+1)=(1-exp(-(i/Nmin)^nPow))*exp(-(i/Nmax10)^2)*exp(-(i/Nmax10/f2tof1)^4);  
%          Filter(i+1) = (1-exp(-(i/Nmin)^4));
    end
    for i=(Nsampl/2) : (Nsampl-1), Filter(i+1)=(1-exp(-((Nsampl-i)/Nmin)^nPow))*exp(-((Nsampl-i)/Nmax10)^2)*exp(-((Nsampl-i)/Nmax10/f2tof1)^4);  end
end

% plot(Filter), xlim([0 200])

% for i=1:1024,    filterF10_v2(i)=(1-exp(-(i/Nmin)^2))*exp(-(i/Nmax10)^2)*exp(-(i/Nmax10/f2tof1)^4);  end
% for i=1025:2048, filterF10_v2(i)=(1-exp(-((2048-i)/Nmin)^2))*exp(-((2048-i)/Nmax10)^2)*exp(-((2048-i)/Nmax10/f2tof1)^4);  end