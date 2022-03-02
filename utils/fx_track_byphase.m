function [lct, pkt, errors_pos, inam_ascan, inph_ascan] = fx_track_byphase(Idx_resin, ori_signal, threshold)
% get the ply thickness
P                           = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply                  = P(1);

% find the peaks in inam, as the front- and back- walls
inam_ascan                  = abs(ori_signal);
MinPeakHeight               = threshold * max(inam_ascan);
MinPeakDistance             = length(inam_ascan) / 2;
[pks_walls, Idx_walls]      = findpeaks(inam_ascan, 1:length(inam_ascan), 'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
front_I                     = round(Idx_walls(1));
back_I                      = round(Idx_walls(end));
inph_ascan                  = angle(ori_signal);
% inph_ascan_unwrap           = unwrap(inph_ascan(front_I:back_I));
% phase_seq                   = inph_ascan_unwrap(1)+pi/2: 2 * pi: inph_ascan_unwrap(end);
% I_phase_track               = interp1(inph_ascan_unwrap, 1:length(inph_ascan_unwrap), phase_seq) + front_I;

% set the region of each layer and find the postion where
% the phase is closest to - pi / 2
pkt                         = NaN(1, length(Idx_resin));
lct                         = NaN(1, length(Idx_resin));
lct(1)                      = front_I;
lct(end)                    = back_I;
pkt(1)                      = pks_walls(1);
pkt(end)                    = pks_walls(end);
phase_0                     = inph_ascan(front_I);

% find the back-wall echo by phase at first
idx_layer_n1                = max(front_I, round(Idx_resin(end) - Idx_oneply / 2));
idx_layer_n2                = min(length(inph_ascan), round(Idx_resin(end) + Idx_oneply / 2));
[pkt(end), lct_back]        = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2) - (phase_0 - pi)));
lct(end)                    = lct_back + idx_layer_n1 - 1;

% % reset the idx_resin by updated back surface !!!! comment out if necessary
% Idx_resin                  = linspace(lct(1), lct(end), length(lct));
% % get the ply thickness again
% P                          = polyfit(1:length(Idx_resin), Idx_resin, 1);
% Idx_oneply                 = P(1);

% % reset Idx_resin by unwrapped phases
% nol                        = length(Idx_resin); % number of layers
% Idx_resin                  = [I_phase_track back_I];
% if length(Idx_resin) > nol
%     Idx_resin = Idx_resin(1:nol);
% end

% % debug
% idx_layer_n1s = ones(1, length(Idx_resin));
% idx_layer_n2s = ones(1, length(Idx_resin));

for n_layer = 2:length(Idx_resin) - 1
    idx_layer_n1                 = round(Idx_resin(n_layer) - Idx_oneply / 2);
    idx_layer_n2                 = round(Idx_resin(n_layer) + Idx_oneply / 2);
%     idx_layer_n1                 = round(lct(n_layer - 1) + Idx_oneply / 2);
%     idx_layer_n2                 = round(lct(n_layer - 1) + Idx_oneply * 3 / 2);
%     % define the true location by one-ply TOF
%     idx_layer_n1                 = round(front_I + Idx_oneply/2 + Idx_oneply*(n_layer-2));
%     idx_layer_n2                 = round(front_I + Idx_oneply/2 + Idx_oneply*(n_layer-1));
    [pkt(n_layer), lct(n_layer)] = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2) - (phase_0-pi/2)));
    lct(n_layer)                 = lct(n_layer) + idx_layer_n1 - 1;
    % check the replicated track points.
    if lct(n_layer)==lct(n_layer-1)
%         disp('replicate');
%         lct(n_layer)             = lct(n_layer) + 1;
        lct(n_layer)             = NaN;
%         lct(n_layer)             = (idx_layer_n1+idx_layer_n2)/2;
    elseif lct(n_layer)>= Idx_resin(end) % exceed the rear surface
        lct(n_layer)             = NaN;
        break;
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

% fillna by average
% lct        = round(fx_inpaint_nans(lct, 1));

% % calculate errors
errors_pos = (lct - Idx_resin)/Idx_oneply;

% find the nearest one to calculate the errors
% errors_pos = NaN(1, length(Idx_resin));
% for idx_l = 1:length(Idx_resin)
%     [~, I_err]        = min(abs(lct - Idx_resin(idx_l)));
%     errors_pos(idx_l) = (lct(I_err) - Idx_resin(idx_l))/Idx_oneply;
% end


% % debug
% figure, plot(inph_ascan);
% hold on;
% scatter(idx_layer_n1s, inph_ascan(idx_layer_n1s), 'rd');
% hold on;
% scatter(idx_layer_n2s, inph_ascan(idx_layer_n2s), 'gs');
% hold on;
% scatter(lct, inph_ascan(lct), 'yo', 'fill');

end

