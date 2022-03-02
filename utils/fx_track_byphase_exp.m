function [lct, pkt, errors_pos, inam_ascan, inph_ascan, phase_0, Idx_resin] = fx_track_byphase_exp(layers_num, ori_signal, threshold)

% find the peaks in inam, as the front- and back- walls
inam_ascan              = abs(ori_signal);
MinPeakHeight           = threshold * max(inam_ascan);
MinPeakDistance         = length(inam_ascan) / 2;
[pks_walls, Idx_walls]  = findpeaks(inam_ascan, 1:length(inam_ascan), 'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
front_I                 = round(Idx_walls(1));
back_I                  = round(Idx_walls(end));
inph_ascan              = angle(ori_signal);

% reset the idx_resin by updated back surface 
Idx_resin               = linspace(front_I, back_I, layers_num);

% set the region of each layer and find the postion where
% the phase is closest to - pi / 2
pkt                     = NaN(1, length(Idx_resin));
lct                     = NaN(1, length(Idx_resin));
lct(1)                  = front_I;
lct(end)                = back_I;
pkt(1)                  = pks_walls(1);
pkt(end)                = pks_walls(end);
phase_0                 = inph_ascan(front_I);

% get the ply thickness again
P                       = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply              = P(1);

% find the back-wall echo by phase at first
idx_layer_n1            = max(front_I, round(Idx_resin(end) - Idx_oneply / 2));
idx_layer_n2            = min(length(inph_ascan), round(Idx_resin(end) + Idx_oneply / 2));
[pkt(end), lct_back]    = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2) - (phase_0 - pi)));
lct(end)                = lct_back + idx_layer_n1 - 1;


% % debug
% idx_layer_n1s = ones(1, length(Idx_resin));
% idx_layer_n2s = ones(1, length(Idx_resin));

for n_layer = 2:length(Idx_resin) - 1
    idx_layer_n1                 = round(Idx_resin(n_layer) - Idx_oneply / 2);
    idx_layer_n2                 = round(Idx_resin(n_layer) + Idx_oneply / 2);
    [pkt(n_layer), lct(n_layer)] = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2) - (phase_0-pi/2)));
    lct(n_layer)                 = lct(n_layer) + idx_layer_n1 - 1;
    % check the replicated track points.
    if lct(n_layer)==lct(n_layer-1)
        lct(n_layer)             = lct(n_layer) + 1;
        if lct(n_layer)==lct(n_layer-1)
            disp('replicate');
        end
    end
%     % debug
%     idx_layer_n1s(n_layer)       = idx_layer_n1;
%     idx_layer_n2s(n_layer)       = idx_layer_n2;
end
% % find the back-wall echo by phase
% idx_layer_n1                     = round(Idx_resin(end) - Idx_oneply / 2);
% idx_layer_n2                     = min(length(inph_ascan), round(Idx_resin(end) + Idx_oneply / 2));
% [pkt(end), lct(end)]             = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2) - (phase_0 - pi)));
% lct(end)                         = lct(end) + idx_layer_n1 - 1;
% calculate errors
errors_pos            = (lct - Idx_resin) / Idx_oneply;
% for i = 1:(length(errors_pos))
%     pos_temp          = min(abs(Idx_resin(i) - lct));
%     errors_pos(i)     = pos_temp / Idx_oneply;
% end

% % cutting the edge of the errors
% errors_pos(errors_pos > 0.1)  = 0.1;
% errors_pos(errors_pos < -0.1) = -0.1;

% % debug
% figure, plot(inph_ascan);
% hold on;
% scatter(idx_layer_n1s, inph_ascan(idx_layer_n1s), 'rd');
% hold on;
% scatter(idx_layer_n2s, inph_ascan(idx_layer_n2s), 'gs');
% hold on;
% scatter(lct, inph_ascan(lct), 'yo', 'fill');

end

