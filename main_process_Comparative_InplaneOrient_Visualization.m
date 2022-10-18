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
%     medi  = median(C_scan_inam_temp(max(wl)/2+2:end-max(wl)/2-1,...
%         max(wl)/2+2:end-max(wl)/2-1),'all');
%     madi = mad(C_scan_inam_temp(max(wl)/2+2:end-max(wl)/2-1,...
%         max(wl)/2+2:end-max(wl)/2-1), 1, 'all');
% 
%     mean_val = mean(C_scan_inam_temp(max(wl)/2+2:end-max(wl)/2-1,...
%         max(wl)/2+2:end-max(wl)/2-1),'all', 'omitnan');
%     std_val = std(C_scan_inam_temp(max(wl)/2+2:end-max(wl)/2-1,...
%         max(wl)/2+2:end-max(wl)/2-1), 1, 'all', 'omitnan');

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


%% read all matrics and plot
close all;
% 
% % Get a list of all files in the folder, and its subfolders
% S = dir(fullfile(FolderName, '*'));
% N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.
% 
% % bar plot elements
% mean_val = nan(4, 4);
% std_val = nan(4, 4);
% ratio = nan(4, 4);
% k = 0;
% 
% for ii = 1:numel(N)
%     disp(ii);
%     T = dir(fullfile(FolderName, N{ii}, '*.csv')); % improve by specifying the file extension.
%     if isempty(T)
%         continue;
%     else
%         k = k + 1;
%         disp(N{ii})
%     end
%     T = readtable(fullfile(T.folder, T.name));
%     mean_val(:, k) = T.Inarea_mean;
%     std_val(:, k) = T.Inarea_std;
%     ratio(:, k) = T.Inarea_ratio;
% end

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
    ratio(:, k) = T.Inarea_ratio;
end

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
end
ylim([-90 -22.5]);
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

xlim([5 30]);
axis off;

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



%% test MumfoldShah 

close all;
clc;

% show C-scan slice
PropertyName = 'img_hil';
glabal_max = 0;
for ratio = [21.5] / 24
    [~, C_scan_inam_temp, ~] = process.define_parallel_inamCscan(ratio, PropertyName);
    % 2D GF-ID
    wavelength  = (4:2:16)*2.5;
    orientation = 1:1:180;
    SFB         = 1; % [0.5 2.5]
    SAR         = 0.5; % [0.23 0.92]
    imagename = 'C_scan_inam_denoise';
    process.(imagename) = C_scan_inam_temp;
    process    = process.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
    % K: controls how much smoothing is applied to the Gabor magnitude responses.
    K = 0.5e0;
    process.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens);
    
    % mumfordshah smooth
    [m, n, k] = size(C_scan_inam_temp); % k = channels, 1 for gray scalar
    % set parameters
    alpha = 5;
    lambda = 0.1;
    eps = 1e-1;
    %
    u_n = fx_mumfordshah(C_scan_inam_temp, alpha, lambda, eps);
    
    % 2D GF-ID
    process.(imagename) = u_n;
    process    = process.compute_logGabor_filter_withoutFig(PropertyName, wavelength, orientation, SFB, SAR, imagename);
    % K: controls how much smoothing is applied to the Gabor magnitude responses.
    K = 1e0;
    process.show_orientation_by_ID_allwl(wavelength, orientation, K, imagename, angle_compens); 
end


