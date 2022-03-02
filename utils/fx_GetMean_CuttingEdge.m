function [mean_errors, std_errors] = fx_GetMean_CuttingEdge(errors_repeat, edge_value)
% calculate the mean value and std
% set the values of the mean value that exceed the edge as the edge_value;

mean_errors = mean(errors_repeat, 1, 'omitnan');
std_errors  = std(errors_repeat, 1, 'omitnan');

mean_errors(mean_errors > edge_value)   = edge_value;
mean_errors(mean_errors < -edge_value) = -edge_value;

end

