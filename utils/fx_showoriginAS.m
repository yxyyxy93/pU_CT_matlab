function  fx_showoriginAS(t, S, fs, fignum)
% t: time domain vector
% S: AS signal vector
% fignum: figure number
% decomposed analytical signal

AS = hilbert(S);

figure(fignum);
subplot(2, 3, 1);
plot(t, S, '-', 'Linewidth', 0.5);
title((' origin signal'), 'fontsize', 16, 'Fontname','times new Roman');
grid on;
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}', 'fontsize', 16);

subplot(2, 3, 2);
plot(t, angle(AS), 'b-');
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}', 'fontsize', 16);

subplot(2, 3, 3);
plot3(t, real(AS), imag(AS), '-', 'Linewidth', 0.5);
grid on;
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}real', 'fontsize', 16);
zlabel('\fontname {times new roman}imag', 'fontsize', 16);
view([90, 0, 0]);

subplot(2, 3, 4);
L = length(S);
Y = fft(S);
P2 = abs(Y/L);
P1 = P2(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);
f = fs*(0:(L/2))/L;
plot(f, P1) 
title('Single-Sided Amplitude Spectrum', 'fontsize', 16, 'Fontname','times new Roman')
xlabel('f (Hz)')
ylabel('|P1(f)|')
hold on;

subplot(2, 3, 5);
plot(t, S);
hold on;
plot(t, abs(AS), 'r-');
xlabel('\fontname {times new roman}t (s)', 'fontsize', 16);
ylabel('\fontname {times new roman}', 'fontsize', 16);
title('real component', 'fontsize', 16, 'Fontname', 'times new Roman');
hold on;

end

