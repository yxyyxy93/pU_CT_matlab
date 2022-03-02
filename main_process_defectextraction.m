% read the file
% used to estimate thickness as well
% This is used for the test the orignal phase-derived interply track and the log-gabor fitelr
% for 1st paper parametric study

clc;
close all;
fclose all;
clear;

%% Define the class
% read the preprocessed .mat
[Filename1, Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename               = strcat(Pathname1, Filename1);
process3               = class_process_woven_RFdata(filename);

% input the settings 
process3.fs = 250e6; % MHz
process3.fx = 1 / 0.2e-3; % unit: 1 / m
process3.fy = 1 / 0.2e-3; % unit: 1 / m

% % read the settings from excel
% [Filename1, Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
% settinig_file          = strcat(Pathname1, Filename1);
% 
% % read the data from the struct, or reset the dataset
% process3               = process3.loadSettings(settinig_file);

process3               = process3.read_origin_data; % read (reset) the dataset

%% ***********
FolderName = "F:\Xiayang\results\Woven_samples\BP1_3_025m\";   % the destination folder
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, 'V313\20220207_BP1_3_025m_V313_25db_PEmat.mat'));

% 20220207_BP1_3_025m_V313_25db_PEmat.mat
% 20220207_BP1_3_025m_V313_25db_DTmat
%% remove DC component and apply hilbert

PropertyName = 'img';
process3     = process3.removeslope(PropertyName);

%% EDA
process3.show_img_shape;
% cut the img, otherwise it could be out of the memory
% window = [101 180 81 160 1 1500]; 
% window = [1 120 151 300 1 1500]; 
window   = [1 400 1 400 200 1200];  
x_step   = 1;
y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);
process3.show_img_shape;
% 3D viewer

volumeViewer(abs(process3.img_hil));

%% show A scan

% define the index to select Ascan 
x           = 45;
y           = 45;
process3.demo_AS_3D_inclinfq(x, y);

%%
% Bscan origin

B_type = 'y';
index = 100;
Bwin  = 1:350;
PropertyName = 'img_hil';
process3.demo_Bscan_inst(B_type, index, Bwin, PropertyName);

%% show C scan
% define the index to select Ascan 
z            = 300;
PropertyName = 'img_hil';
process3.show_Cscan(z, PropertyName);

%% amplitude drop method to determine the defect size

PropertyName = 'img_hil';
filtertype   = {'LP', 'ideal'};
% filtertype   = {'BP', 'ideal'};
% filtertype   = {'nofilter', 'ideal'};

% zrange       = [300 600]; % for DT
zrange       = [900 1200]; % for PE

drop         = -3; % dB
drop         = -6; % dB
process3     = process3.amplitude_drop_method(PropertyName, zrange, drop, filtertype);
 
%% band-pass method to explore the woven structure

filtertype   = {'BP', 'Gaussian'};
% type of the filter: 'ideal', 'Butterworth', or 'Gaussian'

zrange       = [300 600]; % for DT
% zrange       = [900 1200]; % for PE

drop         = -3; % dB
drop         = -6; % dB
process3     = process3.amplitude_drop_method(PropertyName, zrange, drop, filtertype);

%% deconvolution
% read reference signal
[Filename1, Pathname1] = uigetfile({'*.tdms'},  'select the file');   
filename               = strcat(Pathname1, Filename1);
x                      = 10;
y                      = 10;
process3               = process3.read_refer(filename, x, y);

% cut the ref. signal
win                    = [0.1e-6 2e-6]; % unit: s
process3               = process3.cut_reference_signal(win);

% calculate the SNR
win_noise              = [0.1e-6 0.5e-6];
win_signal             = [0.7e-6 1.1e-6]; % unit: s
snr                    = process3.calculate_SNR(win_signal, win_noise);

% align the signal
process3               = process3.align_refer_ascan(x, y); 

%% low-pass filter
bandpassFreq = 30e6;
x = 163;
y = 268;
[ori_signals{1}, t_space, fss] = process3.demo_Ascan(x, y, '2-CEH105-24-p4_50MHz');

process3 = process3.Filter_lowpass(bandpassFreq, 'img_hil');


%% deconvolution
x        = 100;
y        = 100;
q_factor = 1e-2;
ker_wid  = 200;

% process2      = process2.apply_deconvolutions(q_factor, ...
%     q_factor_AR, fft_padding, bw, bwr, k, ker_wid);
process3      = process3.apply_deconvolutions_onlyWiener(q_factor, ker_wid);

% save('process2.mat','process2');

%% surface calculation
% x: x index
% y: y index
% MinPD: MinPeakDistance for findpeaks
% MinPH: MinPeakHeight for findpeaks

PropertyName = 'img_hil_filter';
MinPD        = 25;
MinPH        = 0.005;  % these 2 parameters need to be changed for surface estimation.
alpha        = 4e-3;
% A_ratio      = 0.9;
A_ratio      = 0.95;
% max_len  = 1499;
max_len  = 1250;

% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y
x = 87;
y = 75;

process3.show_Ascan_inam_peaks(x, y, MinPD, MinPH, PropertyName); % x, y
%
process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, x, y)
   
% %
PropertyName = 'img_hil';
PropertyName = 'img_WienerDeconv';
MinPD        = 25;
MinPH        = 0.003;  % these 2 parameters need to be changed for surface estimation.
alpha        = 3e-3;
% A_ratio      = 0.9;
A_ratio      = 0.95;
% max_len  = 1499;
max_len  = 1250;

% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y
x = 163;
y = 268;

process3.show_Ascan_inam_peaks(x, y, MinPD, MinPH, PropertyName); % x, y
%
process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, x, y)
   
%%
% surface calculation
% % the attenuation coefficient for defining the echo from delamiation.
% Ab         = 0.075;
% A0         = 0.83;
% h          = 1156 - 240; % data point 
% z1         = 1000 * 1500;
% z2         = 1588 * 2906;
% [alpha, A] = fx_calcualate_attcoef(Ab, A0, h, z1, z2);
% process3   = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A);
% filename_fig = filename;
% process3.show_surfaces(filename_fig(1:end-5));

% get the surfaces by statistic calculation again
% window       = [1 100 1 100 1 1250];
% [Ab, A0, h]  = process3.statistic_front_backechoes(window);
% z1           = 1000 * 1500;
% z2           = 1588 * 2906;
% 
% [alpha, A]   = fx_calcualate_attcoef(Ab, A0, h, z1, z2);

% alpha        = 5e-3;
% A_ratio      = 0.9;

% A_ratio      = 0.4; % for 0.3 m drop
process3     = process3.find_front_amp_alpha(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
process3.show_surfaces(filename_fig(1:end-5));

% process3     = process3.recover_surface;
% % process3     = process3.smooth_rear_I;
% process3.show_surfaces(filename_fig(1:end-5));

% tackle the problem of walls determination
TOF_walls    = -mean(process3.front_I-process3.rear_I, 'all', 'omitnan');
% process3.rear_I = process3.front_I + mean(TOF_walls);
% process3.show_surfaces;
% % shift the A scan in the time domain
% process3        = process3.shift_A_scan(200);

%% internal damage features
process3.damage_imaging;
%
mid_x     = 178;
mid_y     = 178;
winy      = round(mid_x-23.58e-3*process3.fx/2): round(mid_x+23.58e-3*process3.fx/2);
winx      = round(mid_y-18.67e-3*process3.fy/2): round(mid_y+18.67e-3*process3.fy/2);
% rot_angle = -110;
rot_angle = 20;
process3.show_front_position(winx, winy, rot_angle);

% % 
% centers = [177, 179];
% radii   = 18e-3*process3.fy/2 - 1;
% process3.show_front_position_circle(centers, radii);

%% display A B scans and time gate to determine the back wall

process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, 75, 320);

show_Raw_filtered_3D(obj, xslice, yslice, zslice)
        
process3.demo_Bscan_interply('x', 260, 1:270, 'healthy');

process3.demo_Bscan_interply('x', 160, 1:146, 'delaminated');

volumeViewer(angle(process3.img_hil_filter));

%% monogenic results or filtered

PropertyName = 'img_hil_filter';
process3.find_front_amp_alpha_Ascan(MinPD, MinPH, PropertyName, max_len, alpha, A_ratio, 54, 330);

%% filtered
f0        = 6.5e6;
sigma     = 0.7;

% for debug
% threshold = 0.05;
x         = 54;
y         = 330;
% process3.demo_logGabor_plytrack_inph_v2(x, y, f0, sigma, threshold);
process3.show_hilbert_Ascan(x, y);
% process3.show_logGabor_Scaleogram(3e6:0.1e6:16e6, sigma, x, y);

% 0.676 us - 2.13 us
%  
process3 = process3.Filter_logGabor(f0, sigma, 'img');

%% 3d ply track
% add noise to img.
% process = process.addnoise(25);
% origin without filter
threshold  = 0.01;  % max(inam) * threshold for distinguishing the back surface

% % process  = process.track_interply('img_hil');
% process3   = process3.track_interply_inph(threshold, 'img_hil');
% 
% %
% % process.show_track_interply(xslice, yslice, zslice);
% process3.show_oneinterply(5, 'nofilter', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply(15, 'nofilter', 3000, win_x, win_y, TOF_oneply);
% process3.show_oneinterply(23, 'nofilter', 3000, win_x, win_y, TOF_oneply);
% process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

% % filtered
f0        = 6.5e6;
sigma     = 0.7;

% for debug
% threshold = 0.05;
x         = 100;
y         = 98;
% process3.demo_logGabor_plytrack_inph_v2(x, y, f0, sigma, threshold);
process3.show_hilbert_Ascan(x, y);
process3.show_logGabor_Scaleogram(3e6:0.1e6:16e6, sigma, x, y);

% 0.676 us - 2.13 us
%  
process3 = process3.Filter_logGabor(f0, sigma, 'img');
nol      = 49;

% process3 = process3.track_interply_2ndharmonic('img_hil_filter', nol);
process3 = process3.track_interply('img_hil_filter');

% use 2nd-harmonic for first and the last interplies, 
% use fundamental resonance for another interplies
f0_1     = 12.8e6;
sigma0_1 = 0.8;
f0_2     = 6.3e6;
sigma0_2 = 0.7;
nol      = 26; % make it 25+1 in case that there is redundant accidently-tracked interply
process3 = process3.track_interply_hybrid( ...
    'img', f0_1, sigma0_1, f0_2, sigma0_2, nol);

%%
B_type     = 'y';
index      = 230;
Bwin       = 1:420;
TOF_oneply = TOF_walls / 24; % 24 plies
process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

% 3d slices
xslice = 170 / process3.fx * 1e3;
yslice = 190 / process3.fy * 1e3;
xslice = 175 / process3.fx * 1e3;
yslice = 185 / process3.fy * 1e3;
zslice = [];
process3.show_interply_track_3D_knn(xslice, yslice, zslice);

%% 3D monogenic signal analysis
clc;
D         = process3.img;
[Y, X, Z] = size(D);

cw_xy = 5;
cw_z  = 250e6/15e6;

% woven strucutre
cw_xy = 40;
cw_z  = 250e6/15e6;

% Now use these wavelengths to create a structure containing
% frequency-domain filters to calculate the monogenic signal. Same as for
% the 2D case, except now we need to pass three image dimensions and the
% two sets of wavelengths
filtStruct = createMonogenicFilters3D(Y, X, Z, cw_xy, cw_z, 'lg', 0.55);
% Now we can use this structure to find the monogenic signal for the volume
[m1,m2,m3,m4] = monogenicSignal3D(D, filtStruct);
% Local energy (calculated on a per-scale basis)
LE = localEnergy3D(m1,m2,m3,m4);
% Local amplitude (calculated on a per-scale basis)
LA = sqrt(LE);
% Local phase (calculated on a per-scale basis)
% LP = localPhase3D(m1,m2,m3,m4);
% [FS,FA] = featureSymmetry3D(m1,m2,m3,m4);

clear m1 m2 m3 m4;

% Cscan image
z  = 400;
im = LA(:, :, z);
imagesc(im), axis image; colormap(jet);
h = colorbar;
set(get(h, 'Title'), 'string', '\fontname {times new roman}\fontsize {16} Amp. (arb.)');

% Ascan
x = 45;
y = 45;
figure,plot(squeeze(LA(x, y, :)));

volumeViewer(LA);

%% 3D monogenic signal display
% Display one slice
% 
zslice = 250; % near the middle
figure()
imshow(D(:,:,zslice)), axis image, colormap jet
title('Test Volume Slice')

figure()
imagesc(reshape(LE(:,:,zslice,:),Y,[])), axis image, colormap jet
title('Local Energy Over Scales')

figure()
imagesc(reshape(LP(:,:,zslice,:),Y,[])), axis image, colormap gray
title('Local Phase Over Scales')

figure()
imagesc([FS(:,:,zslice),FA(:,:,zslice)]), axis image, colormap jet
title('3D Feature Symmetry and Asymmetry');

% 
close all;
xslice = 150; % near the middle
figure()
imgshow = squeeze(D(:,xslice,:));
imshow(rot90(imgshow, 3)), axis image, colormap jet
title('Test Volume Slice')

figure()
imgshow = reshape(LE(:,xslice,:,:),Y,[]);
imagesc(rot90(imgshow, 3)), axis image, colormap jet
title('Local Energy Over Scales')

figure()
imgshow = reshape(LP(:,xslice,:,:),Y,[]);
imagesc(rot90(imgshow, 3)), axis image, colormap gray
title('Local Phase Over Scales')

figure()
imgshow = [squeeze(FS(:,xslice,:)),squeeze(FA(:,xslice,:))];
imagesc(rot90(imgshow, 3)), axis image, colormap jet
title('3D Feature Symmetry and Asymmetry');

%% save all figures"z
% "F:\Xiayang\results\Woven_samples\01032021_12h31m02s\f_5_omega07"
% "F:\Xiayang\results\Woven_samples\01032021_12h31m02s\f_10_omega07"
FolderName = "F:\Xiayang\results\Woven_samples\01032021_12h31m02s\f_10_omega07";   % the destination folder
FigList = findobj(allchild(0), 'flat', 'Type', 'figure');
for iFig = 1:length(FigList)
  FigHandle = FigList(iFig);
  FigName   = get(FigHandle, 'Name');
  disp(fullfile(FolderName, FigName));
  set(0, 'CurrentFigure', FigHandle);
  saveas(gcf, strcat(FolderName, '\', FigName), 'epsc');
  saveas(gcf, strcat(FolderName, '\', FigName), 'pdf');
  saveas(gcf, strcat(FolderName, '\', FigName), 'fig');
end

% saveas(gcf,'Barchart','epsc')