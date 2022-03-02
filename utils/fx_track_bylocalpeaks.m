function [lct, pkt, errors_pos, inam] = fx_track_bylocalpeaks(Idx_resin, ori_signal, threshold)
% find the peaks from the original signal by a rule
% get the ply thickness
P             = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply    = P(1);
inam          = abs(hilbert(ori_signal));
% set the region of each layer and find the maximum
pkt           = NaN(1, length(Idx_resin));
lct           = NaN(1, length(Idx_resin));
MinPeakHeight = threshold * max(inam);
[~, locs]     = findpeaks(inam, 'MinPeakHeight', MinPeakHeight);
for n_layer = 1:length(Idx_resin)
    idx_layer_n1                  = round(Idx_resin(n_layer) - Idx_oneply / 2);
    idx_layer_n2                  = round(Idx_resin(n_layer) + Idx_oneply / 2);
    locs_n                        = locs(locs >= idx_layer_n1 & locs < idx_layer_n2);
    if isnan(locs_n)
        lct(n_layer)              = round(idx_layer_n1);
        pkt(n_layer)              = inam(idx_layer_n1);
    elseif ~isnan(locs_n)
        [pkt(n_layer), lct_local] = max(inam(locs_n));
        lct(n_layer)              = round(locs_n(lct_local));
    end
    if isnan(lct(n_layer))
        lct(n_layer)              = round(idx_layer_n1);
    end
end
% calculate errors
errors_pos = (lct - Idx_resin) / Idx_oneply;

end

