% read the file
clc;
close all;
fclose all;
clear;

%%  ************************************ visualization part

close all;
clc;

% show A scan
x   = 50;
y   = 50;
[Ascan, t_space, fss] = process.demo_Ascan(x, y, strcat('A_scan_', 'x'));
xlim([0.5 5.5]);
ylim([-0.4 0.75]);

% zoom in
xlim([4 5]);
ylim([-0.03 0.03]);
xlabel '';
ylabel '';
xticks([])
yticks([])
legend off;

%% 
close all;
clc;

% show C-scan slice
PropertyName = 'img_hil';
glabal_max = 0;
for ratio = [6.5 13.5 17.5 21.5] / 24
    [~, C_scan_inam_temp, ~] = process.define_parallel_inamCscan(ratio, PropertyName);
    x  = (1:size(C_scan_inam_temp, 1))/process.fx*1e3;
    y  = (1:size(C_scan_inam_temp, 2))/process.fy*1e3; 
    cf = figure('Name', strcat('Cscan_Amp_', num2str(ceil(ratio*24))));
    set(cf, 'Position', [0, 0, 500, 400], 'color', 'white');
    pcolor(x, y, C_scan_inam_temp);
    shading flat;
    colormap("jet");
    colorbar;
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
%     clim([0 0.06]);
    if ratio ~= 21.5 / 24
        axis off;
    end
    % save .bmp image to save 
    % glabal_max = max(glabal_max, max(C_scan_inam_temp(:)));
    glabal_max = max(C_scan_inam_temp(:));
    C_scan_inam_temp = round(1 + (C_scan_inam_temp - min(C_scan_inam_temp(:))) ...
        * 255 / (glabal_max - min(C_scan_inam_temp(:))));  % using double
    cmap = gray;
    imwrite(C_scan_inam_temp, cmap, strcat(FolderName, Filename1(1:end-4), ...
        '\Cscan_slice', num2str(round(ratio*24)), '.bmp'))
end

%%
close all;

FolderName = 'C:\Users\xiayang\OneDrive - UGent\matlab_work\results\Different_freq_orientation\';

% show 2D orientation slice
PropertyName = 'Inplane_direction_3D_ID';
% angle_compens = -22; % for 50 MHz
angle_compens = -25;
wl = max(wavelength);

metrics_val = nan(4, 3);
cnt = 0;

for ratio = [6.5 13.5 17.5 21.5] / 24
    [~, C_scan_inam_temp, ~] = process.define_parallel_inamCscan(ratio, PropertyName); % give a modification to not take mean
    x  = (1:size(C_scan_inam_temp, 1))/process.fx*1e3;
    y  = (1:size(C_scan_inam_temp, 2))/process.fy*1e3;
    C_scan_inam_temp = mod(C_scan_inam_temp + angle_compens, 180) - 90;
    % remove the edges
    C_scan_inam_temp(1:max(wl)/2+1, :) = NaN;
    C_scan_inam_temp(:, 1:max(wl)/2+1) = NaN;
    C_scan_inam_temp(end-max(wl)/2:end, :) = NaN;
    C_scan_inam_temp(:, end-max(wl)/2:end) = NaN;
    cf = figure('Name', strcat('Cscan_Orient_', num2str(ceil(ratio*24))));
    set(cf, 'Position', [0, 0, 500, 450], 'color', 'white');
    pcolor(x, y, C_scan_inam_temp);
    shading flat;
    colormap("hsv");
%     colorbar;
    caxis([-90 90]); % below R2022a
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    hold on;
    if ratio ~= 21.5 / 24
        axis off;
    end
    % the gated mean and std
    % the ratio of the leakage area
    bound_l = -68;
    bound_h = -22;
    C_scan_inam_temp_inarea = C_scan_inam_temp(C_scan_inam_temp>bound_l & C_scan_inam_temp<bound_h);
    inarea_ratio = size(C_scan_inam_temp_inarea) / size(C_scan_inam_temp(:));

    inarea_mean = mean(C_scan_inam_temp_inarea, 'omitnan');
    inarea_std  = std(C_scan_inam_temp_inarea, 1, 'omitnan');

    % save the metrics
    cnt = cnt + 1;
    metrics_val(cnt, 1) = inarea_ratio;
    metrics_val(cnt, 2) = inarea_mean;
    metrics_val(cnt, 3) = inarea_std;
    
    % save .bmp image to save 
    C_scan_inam_temp = C_scan_inam_temp(max(wl)/2+1:end-max(wl)/2, ...
        max(wl)/2+1:end-max(wl)/2);
    C_scan_inam_temp = round(C_scan_inam_temp + 90) / 180 * 255;  % using double
    cmap = hsv;
    imwrite(C_scan_inam_temp, cmap, strcat(FolderName, Filename1(1:end-4), ...
        '\ori_slice', num2str(round(ratio*24)), '.bmp'));
end

% names = {'mean' 'std' 'med' 'mad'};
names = {'Inarea_ratio' 'Inarea_mean' 'Inarea_std'};
disp(names);
disp(metrics_val);

% write to csv
T = array2table(metrics_val, "VariableNames", names, "RowNames", {'7', '14', '18', '22'});

if ~exist('SNR','var')
    SNR = 100;
end

writetable(T, strcat(FolderName, Filename1(1:end-4), '\Matrics', ...
    '_', num2str(x_step), '_', num2str(SNR), '_', num2str(process.fs/1e6), '.csv'));

%% 3D display

% show slice
xslice        = 200 / process.fx * 1e3;
yslice        = 200 / process.fy * 1e3;
zslice        = [];
mfsize        = [1 1 1];
angle_compens = -25;
process.show_inplane_direction_3D_ID(xslice, yslice, zslice, mfsize, angle_compens);

axis off;
legend off;
colorbar off;

zlim([1 5.5]);

%% ************** calculate the mean fiber angle and its standard deviation *************
%     % need the reference angle to calcualate the mean and std!
process.statistic_angular_distribution(angle_compens);

ylim([1 5.5]);
% clim([0 0.5]);


%% read all matrics and plot - at all frequencies
close all;

FolderName = ['C:\Users\xiayang\OneDrive - UGent\matlab_work\results\' ...
, 'Different_freq_orientation\dif_freq_matrics'];

S = dir(fullfile(FolderName, '*.csv'));

% bar plot elements
num_fn = numel(S);
mean_val = nan(4, num_fn);
std_val = nan(4, num_fn);
ratio = nan(4, num_fn);
for k = 1:num_fn
    disp(k);
    T = readtable(fullfile(S(k).folder, S(k).name));
    mean_val(:, k) = T.Inarea_mean;
    std_val(:, k) = T.Inarea_std;
    ratio_above(:, k) = T.ratio_above;
    ratio_below(:, k) = T.ratio_below;
end

ratio_a_str = strcat(string(round(ratio_above*100, 1)), '%');
ratio_b_str = strcat(string(round(ratio_below*100, 1)), '%');

% plot setttings
disname = {'50 MHz', '25 MHz', '15 MHz', '5 MHz'};

plotshape = {'s', 'd', 'o', '*', 'h'};
colors = {'red', 'green', 'blue', 'magenta', 'black'};

%  ************* plot the errorbar for plies
cf = figure('Name', '');
set(cf, 'Position', [0, 0, 1200, 600], 'color', 'white');
x = (1:4)* (num_fn + 2);

for k = 1:num_fn
    errorbar(x+k,   mean_val(:, k), std_val(:, k), plotshape{k}, ... 
        'color', colors{k}, 'MarkerSize', 10, 'CapSize', 10, 'display', disname{k}); 
    hold on;
    text(x+k, mean_val(:, k)+std_val(:, k)+5, ratio_a_str(:, k), ...
        'vert','bottom','horiz','center', 'fontname', 'Times new roman', 'fontsize', 12);
    box off;
    text(x+k, mean_val(:, k)-std_val(:, k)-10, ratio_b_str(:, k), ...
        'vert','bottom','horiz','center', 'fontname', 'Times new roman', 'fontsize', 12);
    box off;
end
% ylim([-90 -22.5]);
legend;
set(gca, 'fontname', 'Times new roman');
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);

xticks('');

% % *************  plot the ratio bar ***********
% cf = figure('Name', '');
% set(cf, 'Position', [0, 0, 1200, 300], 'color', 'white');
% ratio_str = strcat(string(round(ratio*100, 1)), '%');
% 
% for k = 1:num_fn
%     h=bar(x+k,  ratio(:, k), 0.12, colors{k});
%     set(h, 'FaceAlpha', 0.2);
%     hold on;
%     text(x+k, ratio(:, k), ratio_str(:, k), 'vert','bottom','horiz','center' ...
%         , 'fontname', 'Times new roman', 'fontsize', 12);
%     box off;
% end
% 
% xlim([5 30]);
% axis off;

%% read all matrics and plot - specify the file names
close all;

foldname = ['C:\Users\xiayang\OneDrive - UGent\matlab_work\results\' ...
'Different_freq_orientation\20220825_2-ceh105-24-p5_v313_23db_PE'];

% filenames = {'Matrics_2_100_250.csv', 
%     'Matrics_4_100_250.csv', 
%     'Matrics_6_100_250.csv', 
%     'Matrics_8_100_250.csv', 
%     'Matrics_10_100_250.csv'};

% filenames = {'Matrics_4_35_250.csv', 
%     'Matrics_4_30_250.csv', 
%     'Matrics_4_25_250.csv', 
%     'Matrics_4_20_250.csv'};

filenames = {
%     'Matrics_4_100_200.csv', 
    'Matrics_4_20_150.csv', 
    'Matrics_4_20_100.csv', 
    'Matrics_4_20_50.csv',
    'Matrics_4_20_30.csv'
    };

% plot setttings
% disname = {'0.1 mm', '0.2 mm', '0.3 mm', '0.4 mm', '0.5 mm'};
% disname = {'35 dB', '30 dB', '25 dB', '20 dB'};
disname = {
%     '200 MS/s', 
    '150 MS/s', '100 MS/s', '50 MS/s', '30 MS/s'};

plotshape = {'s', 'd', 'o', '*', 'h'};
colors = {'red', 'green', 'blue', 'magenta', 'black'};
% bar plot elements
num_fn = numel(filenames);
mean_val = nan(4, num_fn);
std_val = nan(4, num_fn);
ratio = nan(4, num_fn);
k = 0;

for k = 1:num_fn
    disp(k);
    T = readtable(fullfile(foldname, filenames{k}));
    mean_val(:, k) = T.Inarea_mean;
    std_val(:, k) = T.Inarea_std;
    ratio(:, k) = T.Inarea_ratio;
end

%  ************* plot the errorbar for plies
cf = figure('Name', '');
set(cf, 'Position', [0, 0, 1200, 600], 'color', 'white');
x = (1:4)* (num_fn + 2);

for k = 1:num_fn
    errorbar(x+k,   mean_val(:, k), std_val(:, k), plotshape{k}, ... 
        'color', colors{k}, 'MarkerSize', 10, 'CapSize', 10, 'display', disname{k}); 
    hold on;
end
% ylim([-90 -22.5]);
ylim([-70 -30]);
legend;
set(gca, 'fontname', 'Times new roman');
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);

% ************* 
cf = figure('Name', '');
set(cf, 'Position', [0, 0, 1200, 300], 'color', 'white');
ratio_str = strcat(string(round(ratio*100, 1)), '%');

for k = 1:num_fn
    h=bar(x+k,  ratio(:, k), 0.12, colors{k});
    set(h, 'FaceAlpha', 0.2);
    hold on;
    text(x+k, ratio(:, k), ratio_str(:, k), 'vert','bottom','horiz','center' ...
        , 'fontname', 'Times new roman', 'fontsize', 12);
    box off;
end

% xlim([5 35]);
xlim([5 30]);
axis off;

%%    ********************  MumfoldShah results interpretion ************
%     ******************** MumfoldShah results interpretion**************
FolderName = 'C:\Users\xiayang\OneDrive - UGent\matlab_work\results\Different_freq_orientation\';

% test MumfoldShah 
close all;
clc;

%         figure, imagesc(squeeze(u_n(:, :, 500)));
% show C-scan slice
PropertyName = 'img_hil';
glabal_max = 0;
img_3D = abs(process.img_hil);
img_3D = img_3D(:,:,300:round(mean(process.rear_I, 'all')));
[~, ~, lz] = size(img_3D);
% 2D GF-ID
wavelength  = (2:2:16) * 2;
orientation = 1:1:180;
SFB = 1;
SAR = 0.5;
gaborArray = gabor(wavelength, orientation, ...
    'SpatialFrequencyBandwidth', SFB, 'SpatialAspectRatio', SAR);
angle_compens = -25;

% Mumfordshah smoothing
%     figure, imagesc(squeeze(img_3D(:, :, 500)));
if gpuDeviceCount
    img_3D = gpuArray(img_3D);
end

%     set parameters
alpha = 2;
lambda = 0.1;
eps = 1e-1;
u_n = fx_mumfordshah(img_3D, alpha, lambda, eps);
img_3D = gather(img_3D);
u_n = gather(u_n);

for ratio = [21.5] / 24
    %     [~, C_scan_inam_temp, ~] = process.define_parallel_inamCscan(ratio, PropertyName);
    % 2D GF-ID
    C_scan_inam_temp = mean(img_3D(:, :, round(ratio*lz):round(ratio*lz)+2), 3);   
    x  = (1:size(C_scan_inam_temp, 1))/process.fx*1e3;
    y  = (1:size(C_scan_inam_temp, 2))/process.fy*1e3;
    cf = figure('Name', strcat('Cscan_Amp_', num2str(ceil(ratio*24))));
    set(cf, 'Position', [0, 0, 500, 400], 'color', 'white');
    pcolor(x, y, C_scan_inam_temp);
    shading flat;
    colormap gray;   
    colorbar;   
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    saveas(cf, strcat(FolderName, Filename1(1:end-4), ...
        '\Cscan_slice_Ori_', num2str(round(ratio*24)), '.bmp'));
    %
    [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_temp, gaborArray);
    gaborMagnitude      = single(gaborMagnitude);
    % normalized the magnitude.
    wl_Array       = repmat(wavelength', length(gaborArray)/length(wavelength), 1);
    BW             = SFB;
    sigmaX         = wl_Array / pi * sqrt(log(2) / 2) * (2^BW + 1)/(2^BW - 1);
    sigmaY         = sigmaX ./ SAR;
    denominator    = permute(2 * sigmaX .* sigmaY * pi, [2, 3, 1]);
    gaborMagnitude = bsxfun(@rdivide, gaborMagnitude, denominator);
    [~, Idx_1d]    = max(gaborMagnitude, [], 3);
    indAng_2d      = floor(Idx_1d/length(wavelength));
    Ang = orientation(1) + indAng_2d *(orientation(2) - orientation(1));
    Ang = mod(Ang + angle_compens, 180) - 90;  
    cf  = figure('Name', strcat('Cscan_Orient_', num2str(ceil(ratio*24))));
    set(cf, 'Position', [0, 0, 500, 400], 'color', 'white');
    pcolor(x, y, Ang);
    shading flat;
    colormap("hsv");
    colorbar;
    caxis([-90 90]); % below R2022a
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    saveas(cf, strcat(FolderName, Filename1(1:end-4), ...
        '\Inplane_2D_Ori_', num2str(round(ratio*24)), '.bmp'));
    % statistic
    % the ratio of the leakage area
    bound_l = -68;
    bound_h = -22;
    Idof_edge  = -90:90;
    Idof_N = histcounts(Ang, Idof_edge);
    Idof_N = circshift(Idof_N, 45, 2) / sum(Idof_N);
    Idof_edge = Idof_edge - 45;
    Idof_n_inarea = Idof_N(Idof_edge>=bound_l & Idof_edge<=bound_h);
    Idof_edge_inarea = bound_l: bound_h;
    % equation for Gaussian distribution
    gaussEqn = 'a * exp(-1/2*((x-mu)/sigma)^2)';
    % use nlinfit to fit Gaussian using Least Squares
    startPoints = [1 -45 1];
    f1 = fit(Idof_edge_inarea.', Idof_n_inarea.', gaussEqn, 'Start', startPoints ...
        , 'Algorithm', 'Trust-Region', 'MaxIter', 1e2, 'TolFun', 1e-6 ...
        , 'Display', 'iter');
    figure,
    bar(Idof_edge(2:end), Idof_N, 'display', 'angle distribution');
    hold on;
    plot(Idof_edge, f1(Idof_edge), '.-', 'LineWidth', 1, ...
        'display', 'Gaussian fitting');
    legend;
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    saveas(gcf, strcat(FolderName, Filename1(1:end-4), ...
        '\Inplane_distribution_Ori_', num2str(round(ratio*24)), '.bmp'));
    
    %
    C_scan_inam_MS = mean(u_n(:, :, round(ratio*lz):round(ratio*lz)+2), 3);
        cf = figure('Name', strcat('Cscan_Amp_', num2str(ceil(ratio*24))));
    set(cf, 'Position', [0, 0, 500, 400], 'color', 'white');
    pcolor(x, y, C_scan_inam_MS);
    shading flat;
    colormap gray;   
    colorbar;   
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    saveas(cf, strcat(FolderName, Filename1(1:end-4), ...
        '\Cscan_slice_MS_', num2str(round(ratio*24)), '.bmp'));
    %
    [gaborMagnitude, ~] = imgaborfilt(C_scan_inam_MS, gaborArray);
    gaborMagnitude      = single(gaborMagnitude);
    % normalized the magnitude.
    gaborMagnitude = bsxfun(@rdivide, gaborMagnitude, denominator);
    [~, Idx_1d] = max(gaborMagnitude, [], 3);
    indAng_2d   = floor(Idx_1d/length(wavelength));
    Ang_MS = orientation(1) + indAng_2d *(orientation(2) - orientation(1));
    Ang_MS = mod(Ang_MS + angle_compens, 180) - 90;
    cf = figure('Name', strcat('Cscan_Orient_', num2str(ceil(ratio*24))));
    set(cf, 'Position', [0, 0, 500, 400], 'color', 'white');
    pcolor(x, y, Ang_MS);
    shading flat;
    colormap("hsv");
    colorbar;
    caxis([-90 90]); % below R2022a
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    saveas(cf, strcat(FolderName, Filename1(1:end-4), ...
        '\Inplane_2D_MS_', num2str(round(ratio*24)), '.bmp'));
    % statistic
    Idof_edge  = -90:90;
    Idof_N = histcounts(Ang_MS, Idof_edge);
    Idof_N = circshift(Idof_N, 45, 2) / sum(Idof_N);
    Idof_edge = Idof_edge - 45;
    Idof_n_inarea = Idof_N(Idof_edge>=bound_l & Idof_edge<=bound_h);
    Idof_edge_inarea = bound_l: bound_h;
    %
    f1 = fit(Idof_edge_inarea.', Idof_n_inarea.', gaussEqn, 'Start', startPoints ...
        , 'Algorithm', 'Trust-Region', 'MaxIter', 1e2, 'TolFun', 1e-6 ...
        , 'Display', 'iter');
    figure,
    bar(Idof_edge(2:end), Idof_N, 'display', 'angle distribution');
    hold on;
    plot(Idof_edge, f1(Idof_edge), '.-', 'LineWidth', 1, ...
        'display', 'Gaussian fitting');
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    saveas(gcf, strcat(FolderName, Filename1(1:end-4), ...
        '\Inplane_distribution_MS_', num2str(round(ratio*24)), '.bmp'));
end


%% visualize 2 - MumfoldShah result - 3D
close all;

FolderName = 'C:\Users\xiayang\OneDrive - UGent\matlab_work\results\Different_freq_orientation\';

angle_compens = -22; % for 50 MHz
% angle_compens = -25;
xslice        = 100 / process.fx * 1e3;
yslice        = 100 / process.fy * 1e3;
zslice        = [];

inph_ex = Inplane_direction;
inph_ex = mod(inph_ex + angle_compens, 180) - 90;
%             inph_ex = medfilt3(inph_ex, mfsize);
y       = (0: size(inph_ex, 1) - 1) / process.fx * 1e3;
x       = (0: size(inph_ex, 2) - 1) / process.fy * 1e3;
% ******************** show ply tracks
cf = figure('Name', ['3d_orientation_2Dfilter_', 'xslice', num2str(xslice(1))]);
set(cf, 'Position', [0, 0, 600, 300], 'color', 'white');
ax = subplot(1, 1, 1);
x_idx = round(xslice * process.fx / 1e3);
y_idx = round(yslice * process.fy / 1e3);

z         = (0: size(inph_ex, 3) - 1) / process.fs * 1e6; % time axis, unit: us
%                 z         = (0: size(inph_ex, 3) - 1)/ obj.fs * 1e3 * 3000/2; % depth axis
[X, Y, Z] = meshgrid(x, y, z);
% remove the edges
wl = max(wavelength);
inph_ex(1:max(wl)/2+1, :, :) = NaN;
inph_ex(:, 1:max(wl)/2+1, :) = NaN;
inph_ex(end-max(wl)/2:end, :, :) = NaN;
inph_ex(:, end-max(wl)/2:end, :) = NaN;
%
h = slice(ax, X, Y, Z, inph_ex , xslice, yslice, zslice);
hold on;
set(h, 'EdgeColor', 'none');
colormap hsv;
caxis([-90 90]);
axis off;
set(ax, 'ZDir', 'reverse');

saveas(cf, strcat(FolderName, Filename1(1:end-4), ...
    '\tomography_MS', num2str(round(xslice)), '.bmp'));

%% visualize 2 - MumfoldShah result - statistic

[~, ~, lz] = size(inph_ex);
%             Idof_edge  = linspace(-90, 90, length(obj.theta_LG));
Idof_edge  = -90:90;
Idof_Ns    = nan(length(Idof_edge)-1, lz);
for i = 1:lz
    image_orientation = inph_ex(:,:,i);
%     image_orientation = mod(image_orientation+angle_compens, 180) - 90;
    Idof_N            = histcounts(image_orientation, Idof_edge);
    Idof_Ns(:,i)      = log(Idof_N/sum(Idof_N));
end
Idof_Ns   = circshift(Idof_Ns, 22, 1);
Idof_edge = Idof_edge - 22;
%             Idof_edge = Idof_edge + 22 * round(180/length(obj.theta_LG));
figure('Name', ['Angle_distribution_3D']);
imagesc(Idof_edge, (1:lz)/process.fs*1e6, Idof_Ns.');
colormap(jet);
caxis([-12 0]);

axis off;

saveas(gcf, strcat(FolderName, Filename1(1:end-4), ...
    '\statistic_MS', '.bmp'));

%% visualize 2 - MumfoldShah result - 2D
close all;

metrics_val = nan(4, 4);
cnt = 0;

for ratio = [6.5 13.5 17.5 21.5] / 24
    C_scan_inam_temp = inph_ex(:, :, round(ratio*lz));
    % remove the edges
    C_scan_inam_temp(1:max(wl)/2+1, :) = NaN;
    C_scan_inam_temp(:, 1:max(wl)/2+1) = NaN;
    C_scan_inam_temp(end-max(wl)/2:end, :) = NaN;
    C_scan_inam_temp(:, end-max(wl)/2:end) = NaN;
    cf = figure('Name', strcat('Cscan_Orient_', num2str(ceil(ratio*24))));
    set(cf, 'Position', [0, 0, 500, 450], 'color', 'white');
    pcolor(x, y, C_scan_inam_temp);
    shading flat;
    colormap("hsv");
%     colorbar;
    caxis([-90 90]); % below R2022a
    set(gca, 'fontsize', 16);
    set(gca, 'linewidth', 1.5);
    hold on;
    if ratio ~= 21.5 / 24
        axis off;
    end
    %
    % save .bmp image to save 
    C_scan_inam_temp = C_scan_inam_temp(max(wl)/2+1:end-max(wl)/2, ...
        max(wl)/2+1:end-max(wl)/2);
    img_C_scan = round(C_scan_inam_temp + 90) / 180 * 255;  % using double
    cmap = hsv;
    imwrite(img_C_scan, cmap, strcat(FolderName, Filename1(1:end-4), ...
        '\ori_slice_MS', num2str(round(ratio*24)), '.bmp'));

    % the gated mean and std 
    % the ratio of the leakage area
    bound_l = -68;
    bound_h = -22;
    
    % statistic
    Idof_n = exp(Idof_Ns(:, round(ratio*lz)));
    Idof_n_inarea = Idof_n(Idof_edge>=bound_l & Idof_edge<=bound_h);
    Idof_edge_inarea = bound_l: bound_h;
    % equation for Gaussian distribution
    gaussEqn = '(1/sigma/2.5) * exp(-1/2*((x-mu)/sigma)^2)';
    startPoints = [-45 2];
%     gaussEqn = 'a * exp(-1/2*((x-mu)/sigma)^2)';
%     startPoints = [1 -45 2];
    % use nlinfit to fit Gaussian using Least Squares
    f1 = fit(Idof_edge_inarea.', Idof_n_inarea, gaussEqn, 'Start', startPoints ...
        , 'Algorithm', 'Trust-Region', 'MaxIter', 1e2, 'TolFun', 1e-6 ...
        , 'Display', 'iter');
    figure,
    plot(f1, Idof_edge(2:end), Idof_n);
    
    % leakage ratio
    fit_curve = f1(Idof_edge);
    edge_below = Idof_edge<bound_l | (Idof_edge>45 & Idof_edge<135) ...
        & ([Idof_n; 0]>fit_curve).';
    edge_above = Idof_edge>bound_h & Idof_edge<45 & ([Idof_n; 0]>fit_curve).';
    % save the metrics - only the 
    cnt = cnt + 1;
    metrics_val(cnt, 1) = sum(Idof_n(edge_above(1:end-1))); % above leakage
    metrics_val(cnt, 2) = sum(Idof_n(edge_below(2:end))); % below leakage
    metrics_val(cnt, 3) = f1.mu;
    metrics_val(cnt, 4) = f1.sigma;
end

% names = {'mean' 'std' 'med' 'mad'};
names = {'ratio_above' 'ratio_below' 'Inarea_mean' 'Inarea_std'};
disp(names);
disp(metrics_val);

% write to csv
T = array2table(metrics_val, "VariableNames", names, "RowNames", {'7', '14', '18', '22'});

if ~exist('SNR','var')
    SNR = 100;
end

writetable(T, strcat(FolderName, Filename1(1:end-4), '\Matrics_MS', ...
    '_', num2str(x_step), '_', num2str(SNR), '_', num2str(process.fs/1e6), '.csv'));