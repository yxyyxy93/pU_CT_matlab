function [f1, f2] = fx_calculate_bw(f, amp_Vf_r, db)
% ********** find the bandwidth
% f: the frequency domain
% amp_Vf_r: the specturm of the input signal, shoud be a reference pulse.
% db: -3db, -6db, ... the threshold of the bw
% - 10*log10(2) = -3db

[MaxData, MaxI] = max(amp_Vf_r);
ratio = 10^(db / 20);

Threshold=MaxData * ratio;

% 
f1 = interp1(amp_Vf_r(1:MaxI), f(1:MaxI), Threshold);
f2 = interp1(amp_Vf_r(MaxI:end), f(MaxI:end), Threshold);

% plot(f,  amp_Vf_r, 'linewidth', 2);
% hold on;
% line([f1; f1], [0; MaxData], 'linestyle', '--');
% hold on;
% line([f2; f2], [0; MaxData], 'linestyle', '--');
% title("spectrum of intial pulse");
end

