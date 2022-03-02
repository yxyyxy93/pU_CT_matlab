function [lct, pkt, errors_pos, inph_ascan] = fx_track_byphase_matched(I_surface, Idx_resin, ori_signal, phase_0)
% get the ply thickness
% I_surface: the index of the surfaces
P                                = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply                       = P(1);
% find the peaks in inam, as the front- and back- walls
front_I                          = I_surface(1); % take the int for index
back_I                           = I_surface(2);
inph_ascan                       = angle(ori_signal);
% set the region of each layer and find the postion where
% the phase is closest to - pi / 2
pkt                              = NaN(1, length(Idx_resin));
lct                              = NaN(1, length(Idx_resin));
lct(1)                           = front_I;
% % debug
% idx_layer_n1s = ones(1, length(Idx_resin));
% idx_layer_n2s = ones(1, length(Idx_resin));

for n_layer = 1:length(Idx_resin)-1 % including the front and back surfaces
    idx_layer_n1                 = round(Idx_resin(n_layer) - Idx_oneply / 2);
    idx_layer_n2                 = min(length(inph_ascan), round(Idx_resin(n_layer) + Idx_oneply / 2));
    [pkt(n_layer), lct(n_layer)] = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2) - (phase_0-pi/2)));
    lct(n_layer)                 = lct(n_layer) + idx_layer_n1 - 1;
    % debug
%     idx_layer_n1s(n_layer)       = idx_layer_n1;
%     idx_layer_n2s(n_layer)       = idx_layer_n2;
end
% find the back-wall echo by phase
idx_layer_n1                     = round(Idx_resin(end) - Idx_oneply / 2);
idx_layer_n2                     = min(length(inph_ascan), round(Idx_resin(end) + Idx_oneply / 2));
[pkt(end), lct(end)]             = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2) - (phase_0-pi/2)));
lct(end)                         = lct(end) + idx_layer_n1 - 1;

% calculate errors
errors_pos                       = (lct - Idx_resin) / Idx_oneply;

% % debug
% figure, plot(inph_ascan);
% hold on;
% scatter(idx_layer_n1s, inph_ascan(idx_layer_n1s), 'rd');
% hold on;
% scatter(idx_layer_n2s, inph_ascan(idx_layer_n2s), 'gs');
% hold on;
% scatter(lct, inph_ascan(lct), 'yo', 'fill');

end

