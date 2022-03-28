function  fx_showAS(t, S, fs)
% t: time domain vector
% S: AS signal vector
% decomposed analytical signal

S_real = real(S);
S_imag = imag(S);

figure;
subplot(2, 3, 1);
plot(t, S_real, 'k-', 'Linewidth', 2);
title(('origin signal'), 'fontsize', 16, 'Fontname','times new Roman');
grid on;
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}amp. (arb.)', 'fontsize', 16);
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 2);
xlim([t(1) t(end)]);

subplot(2, 3, 2);
L = length(S);
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1: floor(L / 2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = fs * (0:(L / 2)) / L;
plot(f, P1, 'Linewidth', 2);
title('single-sided amplitude spectrum', 'fontsize', 16, 'Fontname','times new Roman');
xlim([1 15e6]);
xlabel('f (Hz)');
ylabel('|P1(f)| (arb.)');
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 2);
hold on;

subplot(2, 3, 3);
plot3(t, S_real, S_imag, '-', 'Linewidth', 2);
title(['analytic-signal '],'fontsize',16,'Fontname','times new Roman');
grid on;
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}real', 'fontsize', 16);
zlabel('\fontname {times new roman}imag', 'fontsize', 16);
view([15, -60, -30]);
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 2);
xlim([t(1) t(end)]);

% subplot(2, 3, 2);
% plot3(t, S_real, S_imag, '-', 'Linewidth', 2);
% title(['analytic-signal (front view)'],'fontsize',16,'Fontname','times new Roman');
% grid on;
% xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
% ylabel('\fontname {times new roman}real', 'fontsize', 16);
% zlabel('\fontname {times new roman}imag', 'fontsize', 16);
% view([90, 0, 0]);
% set(gca, 'fontsize', 16);
% set(gca, 'linewidth', 2);

subplot(2, 3, 4);
plot(t, S_real, 'k-', 'Linewidth', 2);
hold on;
plot(t, abs(S), 'r-', 'Linewidth', 2);
legend({'real component', 'instantaneous amplitude'})
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}amp. (arb.)', 'fontsize', 16);
title('instantaneous amplitude', 'fontsize', 16, 'Fontname', 'times new Roman');
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 2);
hold on;
xlim([t(1) t(end)]);

subplot(2, 3, 5);
plot(t, angle(S), 'b-', 'Linewidth', 2);
title('instantaneous phase', 'fontsize', 16, 'Fontname','times new Roman');
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}phase (rad.)', 'fontsize', 16);
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 2);
xlim([t(1) t(end)]);

subplot(2, 3, 6);
% show inst. freq.
[infq, ~, ~, ~] = fx_ht_inFqPhAm(S_real', fs);
plot(t, infq, 'g-', 'Linewidth', 2);
title(['instantaneous frequency'],'fontsize',16,'Fontname','times new Roman');
grid on;
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}frequency (Hz)', 'fontsize', 16);
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 2);
xlim([t(1) t(end)]);


end

