function errors_pos_repeat = fx_replaceNaN_inCol(errors_pos_repeat, threshold)
% check each column, replace the NaN with the maximum in the column if the
% column contains both NaN and normal values
for idx_col = 1: size(errors_pos_repeat, 2)
    if sum(isnan(errors_pos_repeat(:, idx_col))) < size(errors_pos_repeat, 1) * threshold
        errors_pos_repeat_col                                                   = errors_pos_repeat(:, idx_col);
        errors_pos_repeat_col(isnan(errors_pos_repeat_col)) = max(errors_pos_repeat_col);
        errors_pos_repeat(:, idx_col)                                        = errors_pos_repeat_col;
    else
        errors_pos_repeat(:, idx_col)                                        = NaN;
    end
end

end

