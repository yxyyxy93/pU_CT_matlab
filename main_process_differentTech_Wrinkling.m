% read the file
clc;
close all;
fclose all;
clear;

%%
% 15 MHz Analytic-signal analysis
% Define the class
% read the preprocessed .mat

[Filename1, Pathname1] = uigetfile({'*.mat'}, 'select the file');
filename               = strcat(Pathname1, Filename1);
process3               = class_process_wrinkling_RFdata(filename);

% input the settings 
process3.fs = 250e6; % MHz
process3.fx = 1 / 0.2e-3; % unit: 1 / m
process3.fy = 1 / 0.2e-3; % unit: 1 / m

% % read the settings
% [Filename1, Pathname1] = uigetfile({'*.xlsx'}, 'select the file');
% settinig_file          = strcat(Pathname1, Filename1);
% process3               = process3.loadSettings(settinig_file);

% read the data from the struct, or reset the dataset
process3               = process3.read_origin_data; % read (reset) the dataset

x                      = 10;
y                      = 10;

%% window cutting
% window                = [50 149 40 139 1 1500];

x_step   = 1;
y_step   = 1;
window   = [1 250 1 550 100 1300];
% x_step   = 1;
% y_step   = 1;
process3 = process3.cut_edges(window, x_step, y_step);

process3.show_img_shape;

% shift the A scan in the time domain
process3 = process3.shift_A_scan(0);


%% read reference signal
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
win_signal             = [0.6e-6 1e-6]; % unit: s
snr                    = process3.calculate_SNR(win_signal, win_noise);

% align the signal
process3               = process3.align_refer_ascan(x, y); 

%% deconvolution
x           = 100;
y           = 100;
q_factor    = 1e-2;
q_factor_AR = 1e-2;
fft_padding = 2^11;
f1          = 8e6;
f2          = 19e6;
bw          = [f1 f2];
bwr         = -3:-1:-10;
k           = 45;
lam         = 0.01;
ker_wid     = 200;
Nit         = 45;
DownRate    = 1;
fig_subname = '2-CEH105-24-p4_15MHz';
process3.demo_deconvolutions(x, y, q_factor, q_factor_AR, ...
    fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit, fig_subname);

% demo B scan
B_type      = 'x'; % scaning direction 'y'
index       = 50;
Bwin        = 1:100;
% process = process.deconv_ARextrapolation(Bwin, B_type, index, q_factor_AR, fft_padding, bw, bwr, k, ker_wid, fig_subname);

% process.demo_deconvolutions(x, y, q_factor, q_factor_AR, fft_padding, bw, k, lam, ker_wid, DownRate, Nit, fig_subname);
% 
process3     = process3.apply_deconvolutions_Bscan(Bwin, B_type, index, q_factor, q_factor_AR, fft_padding, bw, bwr, k, lam, ker_wid, DownRate, Nit);

process3.demo_deconvolutions_Bscans(B_type, index, Bwin, fig_subname);

%%
% process2      = process2.apply_deconvolutions(q_factor, ...
%     q_factor_AR, fft_padding, bw, bwr, k, ker_wid);
process3      = process3.apply_deconvolutions_onlyWiener(q_factor, ker_wid);

% save('process2.mat','process2');

%% show A scan
% define the index to select Ascan 
x = 45;
y = 369;
% process3.demo_AS_3D_inclinfq(x, y);

process3.show_hilbert_Ascan(x, y);

%% surface search one signal
global min_pks 
min_pks      = 0.03;
PropertyName = 'img_hil';
% process3.show_Ascan_inam_peaks(172, 150, MinPD, MinPH, PropertyName); % x, y
clc;
close all;

MinPD   = 100;
MinPW   = 13;  % these 2 parameters need to be changed for surface estimation.
% surface calculation
max_len = 800;
alpha   = 7e-3;
A_ratio = 0.4;

x = 63;
y = 234;
process3.find_front_amp_alpha_Ascan(MinPD, MinPW, PropertyName, max_len, alpha, A_ratio, x, y);

%% surface search all signals
% A_ratio      = 0.4; % for 0.3 m drop
process3     = process3.find_front_amp_alpha(MinPD, MinPW, PropertyName, max_len, alpha, A_ratio);
filename_fig = filename;
close all;

process3.show_surfaces(filename_fig(1:end-5));

%% original results display
xslice = 50 / process3.fx * 1e3;
yslice = 50 / process3.fy * 1e3;
zslice = []; % us
process3.show_inaminph_3D(xslice, yslice, zslice);

fx_Scrollable_3d_view(angle(process3.img_hil_filter));
process3.show_hilbert_Ascan(100, 100);

%% internal damage features
process3.damage_imaging;
%
xslice = 67 ;
yslice = 454 ;
process3.show_front_position(xslice, yslice);

%% log-gabor filter
% filtered
f0        = 11e6;
sigma     = 0.7;
process3 = process3.Filter_logGabor(f0, sigma, 'img');

% % low-pass filter
% process3 = process3.Filter_lowpass(10e6, 'img_hil');

fx_Scrollable_3d_view(angle(process3.img_hil_filter));
% fx_Scrollable_3d_view(abs(process3.img_hil));

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

% for debug
% threshold = 0.05;
x         = 100;
y         = 500;
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
f0_1     = 15e6;
sigma0_1 = 0.8;
f0_2     = 7.5e6;
sigma0_2 = 0.75;
nol      = 26; % make it 25+1 in case that there is redundant accidently-tracked interply
process3 = process3.track_interply_hybrid( ...
    'img', f0_1, sigma0_1, f0_2, sigma0_2, nol);

%%
B_type     = 'x';
index      = 300;
Bwin       = 1:100;
TOF_oneply = TOF_walls / 24; % 24 plies
process3.show_B_scan_interply(B_type, index, Bwin, 'img_hil');

% 3d slices
xslice     = 100 / process3.fx * 1e3;
yslice     = 100 / process3.fy * 1e3;
zslice     = [];
process3.show_interply_track_3D_knn(xslice, yslice, zslice);

% process.show_track_interply(xslice, yslice, zslice);
win_x      = 1:80;
win_y      = 1:80;
process3.show_oneinterply(05, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(15, 'logGabor', 3000, win_x, win_y, TOF_oneply);
process3.show_oneinterply(22, 'logGabor', 3000, win_x, win_y, TOF_oneply);
%
std_mean_threshold = [1, 0.22]; % modify here 
figname            = 'logGabor';
process3.thickness_estimation_v1(5.5, 24, win_x, win_y, figname);
% plot with uniformly spaced locaitons
threshold = 0.05;
nol       = 25;
TOF       = process3.track_show_B_scan_interply(...
    B_type, 80, Bwin, 'img_hil_filter', threshold, nol, std_mean_threshold);
% 
% % pos. err. calculation 
% % win_x      = 50:250;
% % win_y      = 50:250;
% layers_num = 25;
% process3   = process3.error_pos_calculate(win_x, win_y, layers_num);

% TOF    = process3.track_show_B_scan_interply(B_type, index, Bwin, 'img_hil_filter', threshold);
save TOF.mat TOF_walls;

%% structrul tensor 
% out-of-plane angles

f0_1     = 10e6;
sigma0_1 = 0.7;
process3 = process3.Filter_logGabor(f0_1, sigma0_1, 'img');

% smoothing scale
ds_rate      = 1;
sigma1       = [3, 3, 3];
% integration scale
sigma2       = [3, 3, 3];
PropertyName = 'img_hil_filter';
process3     = process3.structural_tensor(sigma1, sigma2, PropertyName, ds_rate);

% display the angles 3d slices
medf_kernel  = [1, 1, 1];
xslice       = 60 / process3.fx * 1e3;
yslice       = 65 / process3.fy * 1e3;
zslice       = [];
process3.show_angles_ST_zangle(medf_kernel, xslice, yslice, zslice, ds_rate);
%
B_type     = 'x';
index      = 50;
% Bwin       = 1:750;
Bwin       = 50:800;
process3.show_angles_ST_Bscan(B_type, index, Bwin);
%
B_type     = 'y';
index      = 40;
% Bwin       = 1:200;
Bwin       = 1:97;
process3.show_angles_ST_Bscan(B_type, index, Bwin);

%% video of out-of-plane angles
B_type       = 'x';
Bwin         = 1:850;
propertyname = 'angle_z'; 
process3.makeMovie_angles_ST_Bscan(B_type, Bwin, propertyname);


%% ***********
FolderName = "F:\Xiayang\results\Woven_samples\";   % the destination folder
FolderName = "F:\Xiayang\results\Teijin_Wrinkling\";   % the destination folder
save(strcat(FolderName, Filename1(1:end-4), 'mat'), '-v7.3');

load(strcat(FolderName, '11032021_21h00m37s.mat'));
% '04042021_13h28m29s' 15 MHz
% '11032021_21h00m37s.mat' % 5 MHz

%% **************************************
% one-plane EDA
% 2D demo
% total_plies = 24;
% [angle, process] = process.show_one_Radontransform_vitualprofile(5, 0.5, 'img_hil_filter', center, radii, total_plies);

% % filtered
f0       = 15e6;
sigma    = 0.7; 
process3 = process3.Filter_logGabor(f0, sigma, 'img');

% define the C_scan_inam as well
% PropertyName = 'img_WienerDeconv';
PropertyName = 'img_hil_filter';
% PropertyName = 'img_hil';
ratio        = 8.5/24;
process3     = process3.define_parallel_inamCscan...
    (ratio, PropertyName);
%
imagename = 'C_scan_inam';
% curvelet denoising
sigma     = 5e-4;
process3  = process3.curvelet_denoise_inamCscan(sigma, imagename);

% Radontransform
theta         = 1:1:180;
radiis        = 10;
angle_compens = -7;
process3      = process3.compute_orientation_by_RT_correct(radiis, theta, imagename);
process3.show_orientation_by_ID_RT(radiis, theta, imagename, angle_compens);

% 2d log-Gabor fitler
wavelength  = 2:1:16;
orientation = 0:5:180;
SFB         = 1; % [0.5 2.5]
SAR         = 0.25; % [0.23 0.92]
% imagename = 'C_scan_inam_denoise';
% process3     = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% process3.show_orientation_by_ID_allwl(wavelength, orientation, K);
imagename     = 'C_scan_inam';
angle_compens = 0;
process3      = process3.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
% K: controls how much smoothing is applied to the Gabor magnitude responses.
K = 0.5e0;
process3.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);

%% 'distance to front and rear' in-plane orientation extraction
% RT
% PropertyName     = 'img_WienerDeconv'; % deconvolved signal
PropertyName  = 'img_hil_filter'; % orginal signal
% PropertyName  = 'img_hil';

radius        = 10;
theta         = 1:1:180; % degree
ratio         = 0/24:0.2/24:24/24;
sigma_denoise = 0;
tic;
process3      = process3.extract_local_orientation_RT_3D_parallel(...
    PropertyName, radius, theta, ratio, sigma_denoise);
toc;

% show slice
xslice        = 150 / process3.fx * 1e3;
yslice        = 150 / process3.fy * 1e3;
zslice        = []; % us
angle_compens = 0;
process3.show_inplane_direction_3D(xslice, yslice, zslice, angle_compens);

% 2D log-Gabor filter
wavelength         = 4:2:20;
orientation        = 0:5:180;
ratio              = 0/24:0.2/24:24/24;
K                  = 1e0;
% sigma              = 5e-4;
sigma              = 0;
tic;
process3           = process3.extract_local_orientation_3D_parallel_allwl(...
    PropertyName, wavelength, orientation, ratio, K, sigma);
toc;

% show slice
xslice        = 200 / process3.fx * 1e3;
yslice        = 200 / process3.fy * 1e3;
zslice        = [];
mfsize        = [1 1 1];
angle_compens = 0;
process3.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

% 
B_type = 'x';
index  = 150;
Bwin   = 1:200;
process3.show_inplane_direction_3D_ID_Bscan(B_type, index, Bwin, [3 3 3], angle_compens);
%
B_type = 'y';
index  = 50;
Bwin   = 1:300;
process3.show_inplane_direction_3D_ID_Bscan(B_type, index, Bwin, [3 3 3], angle_compens);

%% ************** calculate the mean fiber angle and its standard deviation *************
% need the reference angle to calcualate the mean and std!
Stacking_sequence = [
    45 0 -45 -90 ...
    45 0 -45 -90 ...
    45 0 -45 -90 ...
    -90 -45 0 45 ...
    -90 -45 0 45 ...
    -90 -45 0 45];
m_fiber_angle_arr   = NaN(1, 24);
std_fiber_angle_arr = NaN(1, 24);
Idof_N_arr          = NaN(24, 180);
PropertyName_IFD    = 'Inplane_direction_3D_ID';
for p = 1:24 % 24 layers for default
    if p ==24
        ratio = (p-0.5) * 1/24;
    else
        ratio = (p-0.5) * 1/24;
    end
    ratio_next = (p-0.3) * 1/24;
    %     ref_angle  = Stacking_sequence(p);
    Idof_N_arr(p,:) = process3.calculate_angle_distrib(PropertyName, PropertyName_IFD, ratio, angle_compens);
    %     std_fiber_angle_arr(p) = std_fiber_angle_arr(p) / sum(Idof_Nr_arr(p,:));
    % normalize
    %     Idof_N_arr(p,:) = (Idof_N_arr(p,:)) / sum(Idof_N_arr(p,:));
    %     % log
    %     Idof_N_arr(p,:) = log(Idof_N_arr(p,:));
    disp([num2str(p) '/' '24']);
    if p==3 || p==11 || p==22
        figure('Name', ['angle_distribution_' num2str(p) '_' PropertyName_IFD]);
        set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
        h1 = bar(-89-22:90-22, Idof_N_arr(p,:)/max(Idof_N_arr(p,:)));
        hold on;
        xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
        ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
        %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
        set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
        set(gca, 'linewidth', 2);
%         set(gca, 'YScale', 'log')
        %         legend([h1 h2], 'Angle distribution', 'Multiple Gaussian fitting', 'Interply track');
        %         ytickformat('percentage');
    end
end
%
maps             = cell(1,24);
for i = 1:24
    maps{i} = ['Ply' num2str(i) ':    ' num2str(Stacking_sequence(i)) '\circ'];
end
fx_creat_chartPlot(maps, Idof_N_arr, Stacking_sequence);

% m_fiber_angle_arr = m_fiber_angle_arr-90-23+angle_compens;
% 
m_fiber_angle_arr_RT   = NaN(1, 24);
std_fiber_angle_arr_RT = NaN(1, 24);
Idof_N_arr             = NaN(24, 180);
PropertyName_IFD       = 'Inplane_direction_3D';
angle_compens          = -7;
for p = 1:24 % 24 layers for default
    if p ==24
        ratio = (p-0.5) * 1/24;
    else
        ratio = (p-0.5) * 1/24;
    end
    ratio_next = (p-0.3) * 1/24;
    ref_angle  = Stacking_sequence(p);
    %     [m_fiber_angle_arr_RT(p), std_fiber_angle_arr_RT(p), yhat, Idof_N_arr(p,:)] = process3.calculate_m_std_fiber_angle( ...
    %         PropertyName, PropertyName_IFD, ratio, ratio_next, angle_compens);
    %     std_fiber_angle_arr_RT(p) = std_fiber_angle_arr_RT(p) / sum(Idof_N_arr(p,:));
    Idof_N_arr(p,:) = process3.calculate_angle_distrib(PropertyName, PropertyName_IFD, ratio, angle_compens);
    %
    %     % normalize
%     Idof_N_arr(p,:) = (Idof_N_arr(p,:)) / sum(Idof_N_arr(p,:));
%     % log
%     Idof_N_arr(p,:) = log(Idof_N_arr(p,:));
    disp([num2str(p) '/' '24']);
    if p==3 || p==11 || p==22
        figure('Name', ['angle_distribution_' num2str(p) '_' PropertyName_IFD]);
        set(gcf, 'Position', [0, 0, 400, 250], 'color', 'white');
        h1 = bar(-89-22:90-22, Idof_N_arr(p,:)/max(Idof_N_arr(p,:)));
        xlabel('\fontname {times new roman} Angle (\circ)', 'fontsize', 16);
        ylabel('\fontname {times new roman} Normalized value', 'fontsize', 16);
        %         ylabel('\fontname {times new roman} Percentage', 'fontsize', 16);
        set(gca, 'Fontname', 'times new Roman', 'FontSize', 16);
        set(gca, 'linewidth', 2);
%         set(gca, 'YScale', 'log')
        %         legend([h1 h2], 'Angle distribution', 'Multiple Gaussian fitting', 'Interply track');
        %         ytickformat('percentage');
    end
end
fx_creat_chartPlot(maps, Idof_N_arr, Stacking_sequence);

% m_fiber_angle_arr_RT = m_fiber_angle_arr_RT-90-23 + angle_compens;
% 
% T = array2table([m_fiber_angle_arr; m_fiber_angle_arr_RT; ...
%     std_fiber_angle_arr; std_fiber_angle_arr_RT]);
% writetable(T,'test.xlsx','Sheet',1);

% %% save all figures 4th paper
% FolderName  = "F:\Xiayang\results\Impact_analysis";   % the destination folder
% FigList     = findobj(allchild(0), 'flat', 'Type', 'figure');
% for iFig = 1:length(FigList)
%   FigHandle = FigList(iFig);
%   FigName   = get(FigHandle, 'Name');
%   disp(fullfile(FolderName, FigName));
%   set(0, 'CurrentFigure', FigHandle);
%   saveas(gcf, strcat(FolderName, '\', FigName), 'epsc');
%   saveas(gcf, strcat(FolderName, '\', FigName), 'pdf');
%   saveas(gcf, strcat(FolderName, '\', FigName), 'fig');
% end

