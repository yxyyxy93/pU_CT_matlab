function [lct_C, inam_ascan, Ascan_phase] = fx_track_bycluster(Ascan_as, threshold, k, Idx_resin)
% track the -pi/2 by distance and clustering
Ascan_phase            = angle(Ascan_as);
% find the peaks in inam, as the front- and back- walls
inam_ascan             = abs(Ascan_as);
MinPeakHeight          = threshold * max(inam_ascan);
MinPeakDistance        = length(inam_ascan) / 2;
[~, Idx_walls]         = findpeaks(inam_ascan, 1:length(inam_ascan), ...
    'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
front_Index            = round(Idx_walls(1));
back_Index             = round(Idx_walls(end));
phase_0                = Ascan_phase(front_Index);
lct                    = find(abs(Ascan_phase-(phase_0-pi/2)) < pi/5);
lct(lct<=front_Index)  = [];
lct(lct>=back_Index)   = [];
if size(lct,1) ==1
    lct = lct';
end
[~, lct_C]             = kmeans(lct, k, 'Distance','sqeuclidean', 'MaxIter', 50, 'Start', Idx_resin');
lct_C                  = round(lct_C);

end

