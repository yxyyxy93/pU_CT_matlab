function  fx_showAS(t, S, fs)
% t: time domain vector
% S: AS signal vector
% decomposed analytical signal

S_real = real(S);
S_imag = imag(S);

figure;
subplot(2, 3, 1);
plot(t, S_real, '-', 'Linewidth', 1);
title((' origin signal'), 'fontsize', 16, 'Fontname','times new Roman');
grid on;
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}', 'fontsize', 16);

subplot(2, 3, 2);
plot(t, angle(S), 'b-', 'Linewidth', 1);
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}', 'fontsize', 16);

subplot(2, 3, 3);
plot3(t, S_real, S_imag, '-', 'Linewidth', 1);
grid on;
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}real', 'fontsize', 16);
zlabel('\fontname {times new roman}imag', 'fontsize', 16);
view([90, 0, 0]);

subplot(2, 3, 4);
L = length(S);
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1: floor(L / 2+1));
P1(2:end-1) = 2*P1(2:end-1);
f = fs * (0:(L / 2)) / L;
plot(f, P1, 'Linewidth', 1) 
title('Single-Sided Amplitude Spectrum', 'fontsize', 16, 'Fontname','times new Roman')
xlabel('f (Hz)')
ylabel('|P1(f)|')
hold on;

subplot(2, 3, 5);
plot(t, S_real, 'Linewidth', 1);
hold on;
plot(t, abs(S), 'r-', 'Linewidth', 1);
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}', 'fontsize', 16);
title('real component', 'fontsize', 16, 'Fontname', 'times new Roman');
hold on;

subplot(2, 3, 6);
plot3(t, S_real, S_imag, '-', 'Linewidth', 1);
title(['analytical signal '],'fontsize',16,'Fontname','times new Roman');
grid on;
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}real', 'fontsize', 16);
zlabel('\fontname {times new roman}imag', 'fontsize', 16);
view([15, -60, -30]);

% show inst. freq.
[infq, ~, ~, ~] = fx_ht_inFqPhAm(S_real', fs);
figure,
plot(t, infq, '-', 'Linewidth', 1);

end

