function fx_demo_best_results(Origin_errors_50, WDARextrap_error_pos, AS_LG_errors_pos, edge_value, edge_value_std)
% demonstrate the best result of each of the previous techniques
% edge_value, edge_value_std: the edge values for mean and std
% different markers
Markers_l        = {'-rd' '-go' '-bs' '-ch' '-m*' '-y+' '-kv', '-cx'};
Markers          = Markers_l;

% % substract
% Origin_errors_50     = Origin_errors_50(:, :, 2:end) - Origin_errors_50(:, :, 1:end-1);
% WDARextrap_error_pos = WDARextrap_error_pos(:, :, 2:end) - WDARextrap_error_pos(:, :, 1:end-1);
% AS_LG_errors_pos     = AS_LG_errors_pos(:, :, 2:end) - AS_LG_errors_pos(:, :, 1:end-1);

% % substract front
% Origin_errors_50     = Origin_errors_50(:, :, 2:end) - Origin_errors_50(:, :, 1);
% WDARextrap_error_pos = WDARextrap_error_pos(:, :, 2:end) - WDARextrap_error_pos(:, :, 1);
% AS_LG_errors_pos     = AS_LG_errors_pos(:, :, 2:end) - AS_LG_errors_pos(:, :, 1);
% 
% substract average 
Origin_errors_50     = Origin_errors_50(:, :, 1:end) - mean(Origin_errors_50(:, :, 1:end), 3);
WDARextrap_error_pos = WDARextrap_error_pos(:, :, 1:end) - mean(WDARextrap_error_pos(:, :, 1:end), 3);
AS_LG_errors_pos     = AS_LG_errors_pos(:, :, 1:end) - mean(AS_LG_errors_pos(:, :, 1:end), 3);

% cutting the edge of the errors
Origin_std_50    = squeeze(std(Origin_errors_50, 0, [1 2], 'omitnan'));
WDARextrap_std   = squeeze(std(WDARextrap_error_pos, 0, [1 2], 'omitnan'));
AS_LG_std        = squeeze(std(AS_LG_errors_pos, 0, [1 2], 'omitnan'));
Origin_errors_50 = squeeze(mean(abs(Origin_errors_50), [1 2], 'omitnan'));
WDARextrap_error = squeeze(mean(abs(WDARextrap_error_pos), [1 2], 'omitnan'));
AS_LG_errors     = squeeze(mean(abs(AS_LG_errors_pos), [1 2], 'omitnan'));
%
Origin_errors_50(Origin_errors_50>edge_value)    = edge_value;
Origin_errors_50(Origin_errors_50<-edge_value)   = -edge_value;
WDARextrap_error(WDARextrap_error>edge_value)    = edge_value;
WDARextrap_error(WDARextrap_error<-edge_value)   = -edge_value;
AS_LG_errors(AS_LG_errors>edge_value)            = edge_value;
AS_LG_errors(AS_LG_errors<-edge_value)           = -edge_value;
% plot the errors;
cf2                           = figure('Name', ['Trackerrors_compare_', 'exp']);
set(cf2, 'Position', [0, 0, 800, 1000], 'color', 'white');
ca = subplot(2, 1, 1);
plot(ca, Origin_errors_50, Markers{1}, ...
    'MarkerSize', 8, 'linewidth', 2, 'DisplayName', '50MHz & envelope');
hold on;
plot(ca, WDARextrap_error, Markers{2}, ...
    'MarkerSize', 8, 'linewidth', 2, 'DisplayName', '15MHz & Wiener deconvolution with spectrum extrapolation');
hold on;
plot(ca, AS_LG_errors, Markers{3}, ...
    'MarkerSize', 8, 'linewidth', 2, 'DisplayName', '5 MHz & Analytical signal with log-Gabor filter');
hold on;
xlabel(ca, '\fontname {times new roman}\fontsize {16} Time (\mus)');
ylabel(ca, '\fontname {times new roman} Mean of absolute errors', 'fontsize', 16);
% ytickformat('percentage');
% ylim([0 edge_value*100]);
set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
set(ca, 'linewidth', 2);
legend show;
legend('location', 'best');
legend('fontsize', 12);
grid on;
% cut the edge by edge value
Origin_std_50(Origin_std_50>edge_value_std)       = edge_value_std;
WDARextrap_std(WDARextrap_std>edge_value_std)     = edge_value_std;
AS_LG_std(AS_LG_std>edge_value_std)               = edge_value_std;
%
ca = subplot(2, 1, 2);
plot(ca, Origin_std_50, Markers_l{1}, ...
    'MarkerSize', 8, 'linewidth', 2, 'DisplayName', '50MHz & envelope');
hold on;
plot(ca, WDARextrap_std, Markers_l{2}, ...
    'MarkerSize', 8, 'linewidth', 2, 'DisplayName', '15MHz & Wiener deconvolution with spectrum extrapolation');
hold on;
plot(ca, AS_LG_std, Markers_l{3}, ...
    'MarkerSize', 8, 'linewidth', 2, 'DisplayName', '5 MHz & Analytical signal with log-Gabor filter');
hold on;
xlabel(ca, '\fontname {times new roman}\fontsize {16} Time (\mus)');
ylabel(ca, '\fontname {times new roman} Std of errors', 'fontsize', 16);
% ytickformat('percentage');
% ylim([0 edge_value_std*100]);
set(ca, 'Fontname', 'times new Roman', 'Fontsize', 16);
set(ca, 'linewidth', 2);
legend show;
legend('location', 'best');
legend('fontsize', 12);
grid on;

end