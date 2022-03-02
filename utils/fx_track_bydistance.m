function [lct_C, inam_ascan, Ascan_phase, errors_pos] = fx_track_bydistance(Ascan_as, threshold, Idx_resin)
% get the ply thickness
P                      = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply             = P(1);
% track the -pi/2 by distance
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
% find the back-wall echo by phase at first
idx_layer_n1           = round(Idx_resin(end) - Idx_oneply / 2);
idx_layer_n2           = min(length(Ascan_phase), round(Idx_resin(end) + Idx_oneply / 2));
[~, lct_back]          = min(abs(Ascan_phase(idx_layer_n1: idx_layer_n2) - (phase_0 - pi)));
lct_back               = lct_back + idx_layer_n1 - 1;
% reset the idx_resin by updated back surface
Idx_resin              = linspace(front_Index, lct_back, length(Idx_resin));
% get the ply thickness again
P                      = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply             = P(1);

lct                    = find(abs(Ascan_phase-(phase_0-pi/2)) < pi/5);
lct(lct<=front_Index)  = [];
lct(lct>=back_Index)   = [];

lct_C_temp = lct(1);
lct_C      = [];
for i = 2:length(lct)
    if lct(i) - lct_C_temp(1) <= Idx_oneply / 2
        lct_C_temp = [lct_C_temp lct(i)];
    else
        lct_C = [lct_C mean(lct_C_temp)];
        lct_C_temp = lct(i);
    end
    if i==length(lct)
        lct_C = [lct_C mean(lct_C_temp)];
    end
end

lct_C_all             = [front_Index lct_C lct_back];
errors_pos            = NaN(1, length(lct_C_all));
for i = 1:(length(errors_pos))
    pos_temp          = min(abs(Idx_resin(i) - lct_C_all));
    errors_pos(i) = pos_temp / Idx_oneply;
end

lct_C                 = round(lct_C);

    
% debug
figure, plot(Ascan_phase);
hold on; 
scatter(lct, Ascan_phase(lct));
figure, plot(Ascan_phase)
hold on; 
scatter(lct_C, Ascan_phase(lct_C));
figure, plot(Ascan_phase)
hold on; 
scatter(round(Idx_resin), Ascan_phase(round(Idx_resin)));


end

