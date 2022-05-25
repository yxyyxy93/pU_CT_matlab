% read the file

clc;
close all;
fclose all;
clear;

%% ***********
FolderName = "C:\Users\xiayang\OneDrive - UGent\matlab_work\results\Impact_damaged_sample\";   % the destination folder

load(strcat(FolderName, '2-CEH105-24-p5-6-030m_GEH5M.mat'));

%% check A-scans - 2 times fft
% close all;

% x = 155;
% y = 75;
x = 55;
y = 75;

% temp     = img_noiseadded(x, y, :);
temp     = temp_img(x, y, :);
temp     = squeeze(temp);
Fs       = 250e6; % MHz
% 
L = length(temp);
L = 2^nextpow2(L);
figure,
t_space = 1/Fs:1/Fs:L/Fs;
plot(t_space, [temp; zeros(L-length(temp), 1)], 'LineWidth', 2);
hold on;
plot(t_space, [abs(hilbert(temp)); zeros(L-length(temp), 1)], 'LineWidth', 2);


temp_fft = fft(temp, L)/L;
temp_fft = temp_fft(1:L/2+1);

f_space  = 0:Fs/L:Fs/2;

temp_fft_abs = abs(temp_fft);
figure,
plot(f_space, temp_fft_abs, 'LineWidth', 2);
% figure,
% plot(diff(temp_fft_abs), 'LineWidth', 2);

temp_fft2_signal = fft(temp_fft_abs);
t_space_2        = 0:2/Fs:L/Fs/2;
temp_fft2_signal = temp_fft2_signal(1:length(t_space_2));
figure,
plot(t_space_2, abs(temp_fft2_signal), 'LineWidth', 2);

%% add noise
close all;
dBW      = -25;
[lx, ly, lz] = size(temp_img);
img_noiseadded = zeros(lx, ly, lz);
for i = 1:lx
    for j = 1:ly
        ascan                   = squeeze(temp_img(i, j, :));
        img_noiseadded(i, j, :) = ascan + wgn(lz, 1, dBW);
    end
    clc;
    disp([num2str(i) '/' num2str(lx)]);
end

%% downsample dataset
close all;
rate           = 250/10;
[lx, ly, lz]   = size(temp_img);
img_noiseadded = zeros(lx, ly, round(lz/rate));
for i = 1:lx
    for j = 1:ly
        ascan = downsample(temp_img(i, j, :), round(rate));
        img_noiseadded(i, j, 1:length(ascan)) = ascan;
    end
    clc;
    disp([num2str(i) '/' num2str(lx)]);
end

