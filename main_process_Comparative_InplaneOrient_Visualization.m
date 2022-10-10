% read the file
clc;
close all;
fclose all;
clear;

filevector = {'F:\Xiayang\2-ceh105-24-p5\v390\20220705_2-ceh105-24-p5_v390_25db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\v324\20220825_2-ceh105-24-p5_v324_22db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\v313\20220825_2-ceh105-24-p5_v313_23db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\v320\20220825_2-ceh105-24-p5_v320_17db_PE.mat', ...
    ... 'F:\Xiayang\2-ceh105-24-p5\v306\20220705_2-ceh105-24-p5_v306_21db_PE.mat', ...
    'F:\Xiayang\2-ceh105-24-p5\geh5m\20220825_2-ceh105-24-p5_h5m_8db_PE.mat'};

fsvector      = [250 250 250 250 250 250] * 1e6;
windowsvector = {...
    [234 633 126 525 1 1400], ...
    [56 455 87 486 1 1500], ...
    [56 455 87 486 1 1500], ...
    [56 455 87 486 1 1400], ...
    ... [22 421 95 494 1 1400], ...
    [33 432 81 480 1 1400]};
process.show_surfaces('surface');

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
    clim([0 0.06]);
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
% show 2D orientation slice
PropertyName = 'Inplane_direction_3D_ID';
% angle_compens = -21; % for 50 MHz
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
    clim([-90 90]);
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
    inarea_ratio = size(C_scan_inam_temp_inarea) / size(C_scan_inam_temp(~isnan(C_scan_inam_temp)));

    inarea_mean = mean(C_scan_inam_temp_inarea, 'omitnan');
    inarea_std  = std(C_scan_inam_temp_inarea, 'omitnan');

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

writetable(T, strcat(FolderName, Filename1(1:end-4), '\Matrics', ...
    '_', num2str(x_step), '_', num2str(SNR), '.csv'));

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

% Get a list of all files in the folder, and its subfolders
S = dir(fullfile(FolderName, '*'));
N = setdiff({S([S.isdir]).name},{'.','..'}); % list of subfolders of D.

% bar plot elements
mean_val = nan(4, 4);
std_val = nan(4, 4);
ratio = nan(4, 4);
k = 0;

for ii = 1:numel(N)
    disp(ii);
    T = dir(fullfile(FolderName, N{ii}, '*.csv')); % improve by specifying the file extension.
    if isempty(T)
        continue;
    else
        k = k + 1;
        disp(N{ii})
    end
    T = readtable(fullfile(T.folder, T.name));
    mean_val(:, k) = T.Inarea_mean;
    std_val(:, k) = T.Inarea_std;
    ratio(:, k) = T.Inarea_ratio;
end

%  ************* plot the errorbar for plies
cf = figure('Name', '');
set(cf, 'Position', [0, 0, 1200, 600], 'color', 'white');
x = [1, 2, 3, 4]*6;
errorbar(x,   mean_val(:, 1), std_val(:, 1), ... 
    's', 'MarkerSize', 10, 'CapSize', 10, 'color', 'red', 'display', '50 MHz'); 
hold on;
errorbar(x+1, mean_val(:, 4), std_val(:, 4), ...
    'd', 'MarkerSize',10, 'CapSize', 10, 'color', 'green', 'display', '25 MHz'); 
hold on;
errorbar(x+2, mean_val(:, 3), std_val(:, 3), ...
    'o', 'MarkerSize',10, 'CapSize', 10, 'color', 'blue', 'display', '15 MHz'); 
hold on;
errorbar(x+3, mean_val(:, 2), std_val(:, 2), ...
    '*', 'MarkerSize',10, 'CapSize', 10, 'color', 'magenta', 'display', '5 MHz'); 
hold on;
ylim([-90 -22.5]);
legend;
set(gca, 'fontname', 'Times new roman');
set(gca, 'fontsize', 16);
set(gca, 'linewidth', 1.5);

% ************* 
cf = figure('Name', '');
set(cf, 'Position', [0, 0, 1200, 300], 'color', 'white');
ratio_str = strcat(string(round(ratio*100, 1)), '%');
% h = bar([x', x'+1, x'+2, x'+3],  ratio, 0.8); 
% set(h, 'FaceAlpha', 0.3); 
h=bar(x,  ratio(:, 1), 0.12, 'red');
set(h, 'FaceAlpha', 0.2); 
hold on;
text(x, ratio(:, 1), ratio_str(:, 1), 'vert','bottom','horiz','center' ...
    , 'fontname', 'Times new roman', 'fontsize', 12); 
box off;

h=bar(x+1,  ratio(:, 4), 0.12, 'green');
set(h, 'FaceAlpha', 0.2); 
hold on;
text(x+1, ratio(:, 4), ratio_str(:, 4),'vert','bottom','horiz','center' ...
    , 'fontname', 'Times new roman', 'fontsize', 12); 
box off;

h=bar(x+2,  ratio(:, 3), 0.12, 'blue'); 
set(h, 'FaceAlpha', 0.2); 
hold on;
text(x+2, ratio(:, 3), ratio_str(:, 3), 'vert','bottom','horiz','center' ... 
    , 'fontname', 'Times new roman', 'fontsize', 12); 
box off;

h=bar(x+3,  ratio(:, 2), 0.12, 'magenta'); 
set(h, 'FaceAlpha', 0.2); 
hold on;
text(x+3, ratio(:, 2), ratio_str(:, 2),'vert','bottom','horiz','center' ...
    , 'fontname', 'Times new roman', 'fontsize', 12); 
box off;
xlim([5 30]);
axis off;


%% read all matrics and plot - specify the file names
close all;

foldname = ['F:\Xiayang\results\DifferentTech_fiberorientation\' ...
    'comparison\20220825_2-ceh105-24-p5_v313_23db_PE'];

% filenames = {'Matrics_2.csv', 
%     'Matrics_4.csv', 
%     'Matrics_6.csv', 
%     'Matrics_8.csv', 
%     'Matrics_10.csv'};

filenames = {'Matrics_4_35.csv', 
    'Matrics_4_30.csv', 
    'Matrics_4_25.csv', 
    'Matrics_4_20.csv'};

% plot setttings
% disname = {'0.1 mm', '0.2 mm', '0.3 mm', '0.4 mm', '0.5 mm'};
disname = {'35 dB', '30 dB', '25 dB', '20 dB'};
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
