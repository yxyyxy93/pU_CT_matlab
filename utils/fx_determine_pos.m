function [determined_pos, error_pos] = fx_determine_pos(t_space, ori_signal, threshold, MinPeakDistance, Idx_resin, TOF_oneply)
% use threshold to identify the piles from the ultrassonic response
% calculate the track errors 

determined_pos = NaN(1, length(Idx_resin));
error_pos          = NaN(1, length(Idx_resin));
MinPeakHeight = max(ori_signal) * threshold;
[~, lct]              = findpeaks(abs(ori_signal), t_space, 'MinPeakHeight', MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
% at each theoretical positions of the resin layers, find the nearest locations as the determined positions
for idx = 1:length(Idx_resin)
    [~, I_pos]= min( abs(lct - t_space(Idx_resin(idx))) );
    % dissmiss the duplicated postions.
    if ismember(lct(I_pos), determined_pos)
        continue;
    end
    determined_pos(idx) = lct(I_pos);
    error_pos( idx)  = lct(I_pos) - t_space(Idx_resin(idx));
    % if the determined_pos is more than half of the ply length, dismiss it
    if abs(error_pos(idx)) >= TOF_oneply / 2
        determined_pos(1, idx) = NaN;
        error_pos(1, idx)  = NaN;
    end
end
determined_pos = determined_pos(isfinite(determined_pos)); % dismiss the NaN

end

