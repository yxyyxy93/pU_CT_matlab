function [lct, pkt, errors_pos, inam_ascan, inph_ascan] = fx_track_byphase_2pi_increment(Idx_resin, ori_signal, threshold)
% get the ply thickness
P                 = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply        = P(1);

% find the peaks in inam, as the front- and back- walls
inam_ascan        = abs(ori_signal);
% MinPeakHeight     = threshold * max(inam_ascan);
% MinPeakDistance   = length(inam_ascan) / 3;
% [pks_walls, ~]    = findpeaks(inam_ascan, 1:length(inam_ascan), 'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
front_I           = round(Idx_resin(1));
back_I            = round(Idx_resin(end));
% %do a parabolic estimate
% [front_I, ~, ~]        = crit_interp_p(inam_ascan(front_I+(-1:1)), front_I+(-1:1));
% [back_I, ~, ~]         = crit_interp_p(inam_ascan(back_I+(-1:1)), back_I+(-1:1));

% define the 2\pi phase cycles
inph_ascan        = angle(ori_signal);
tol               = pi;
inph_ascan_unwrap = unwrap(inph_ascan, tol);
if round(back_I) >0 && ...
        round(front_I) > 0 && ...
        inph_ascan_unwrap(round(back_I)) - inph_ascan_unwrap(round(front_I))+pi >= 2*pi
    phase_seq     = inph_ascan_unwrap(round(front_I))+pi: 2*pi: inph_ascan_unwrap(round(back_I))-pi; % back_inph - pi!
    % remove duplicates using 'unique' function.
    [~, ind]      = unique(inph_ascan_unwrap); % ind = index of first occurrence of a repeated value
    v             = 1:length(inph_ascan_unwrap);
    I_phase_track = interp1(inph_ascan_unwrap(ind), v(ind), phase_seq);
else
    I_phase_track = [];
end

% % for debug
% figure, ca = subplot(3,1,1);
% plot(inam_ascan);
% hold on;
% scatter([front_I back_I], inam_ascan(round([front_I back_I])));
% scatter(I_phase_track, inam_ascan(round(I_phase_track)));
% ca = subplot(3,1,2);
% plot(inph_ascan);
% hold on;
% scatter([front_I back_I], inph_ascan(round([front_I back_I])));
% scatter(I_phase_track, inph_ascan(round(I_phase_track)));
% ca = subplot(3,1,3);
% plot(inph_ascan_unwrap);
% hold on;
% scatter([front_I back_I], inph_ascan_unwrap(round([front_I back_I])));
% scatter(I_phase_track, inph_ascan_unwrap(round(I_phase_track)));

% set the region of each layer and find the postion where
% the phase is closest to - pi / 2
pkt                  = NaN(1, length(Idx_resin));
lct                  = NaN(1, length(Idx_resin));
lct(1)               = front_I;
lct(end)             = back_I;

% % interpolation
phase_0              = interp1(round(front_I)+(-1:1), inph_ascan(round(front_I)+(-1:1)), front_I);
% phase_0                 = inph_ascan(front_I);

pkt(1)               = phase_0; % pkt represents the phase

% define phase_0 - pi/2
if phase_0-pi/2 >= -pi
    phase_pi2 = phase_0 - pi/2;
else
    phase_pi2 = phase_0 + pi*3/2;
end

%%
% find the back-wall echo by phase at first
idx_layer_n1         = max(front_I, round(back_I-Idx_oneply/2));
idx_layer_n2         = min(length(inph_ascan), round(back_I+Idx_oneply/2));
if phase_0>0
    [pkt(end), lct_back] = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2)-(phase_0-pi)));
    pkt(end)             = phase_0-pi;
else
    [pkt(end), lct_back] = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2)-(phase_0+pi)));
    pkt(end)             = phase_0+pi;
end
lct(end)             = lct_back+idx_layer_n1-1;

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

%%
% % debug
% idx_layer_n1s = ones(1, length(Idx_resin));
% idx_layer_n2s = ones(1, length(Idx_resin));

% % use for increment of the locations
% pre_loc                  = front_I;
% increment_gap            = 10; % add 10 points gap
   
for n_layer = 2:length(I_phase_track)+1
    %     idx_layer_n1                 = round(Idx_resin(n_layer) - Idx_oneply / 2);
    %     idx_layer_n2                 = round(Idx_resin(n_layer) + Idx_oneply / 2);
    % reaches the true number of the layers, which means this layer should
    % be the back layer
    if n_layer >= length(Idx_resin)
        break;
    end
    idx_layer_n1 = round(I_phase_track(n_layer-1));
    if n_layer==length(I_phase_track)+1
        idx_layer_n2 = round(back_I);
    else
        idx_layer_n2 = round(I_phase_track(n_layer));
    end
    % find one nearest point and apply interpolation
    inph_ascan_oneply         = inph_ascan(idx_layer_n1: idx_layer_n2);
    %in order to only find positive numbers 
    inph_check                = inph_ascan_oneply - phase_pi2;
    inph_check(inph_check<=0) = nan; %replace negative numbers and the zero with nan
    [~, I_minpos]             = min(inph_check); %find values
    %in order to only find negative numbers 
    inph_check                = inph_ascan_oneply - phase_pi2;
    inph_check(inph_check>0)  = nan; %replace negative numbers and the zero with nan
    [~, I_maxneg]             = max(inph_check); %find values
    nearest_ps                = [I_maxneg I_minpos];
    %     % check whether it exceeds the bound
    %     left_b               = max(1, nearest_p-1);
    %     right_b              = min(length(inph_ascan_oneply), nearest_p+1);
    %     nearest_ps           = left_b:1:right_b;
    x                         = inph_ascan_oneply(nearest_ps);
    [~, ind]                  = unique(x); % ind = index of first occurrence of a repeated value
    if isempty(ind)
        lct(n_layer) = NaN;
        pkt(n_layer) = NaN;
    elseif length(ind)==1 || ...
            isempty(find(x>phase_pi2, 1)) || ...
            isempty(find(x<=phase_pi2, 1))
        lct(n_layer) = round(mean(nearest_ps(ind)))+idx_layer_n1-1;
        pkt(n_layer) = phase_pi2;
    else
        vq1                  = interp1(x(ind), nearest_ps(ind), phase_pi2, 'linear');
        lct(n_layer)         = vq1+idx_layer_n1-1;
        pkt(n_layer)         = phase_pi2;
    end
%     % debug
%     idx_layer_n1s(n_layer-1)       = idx_layer_n1;
%     idx_layer_n2s(n_layer-1)       = idx_layer_n2;
end

% % find the back-wall echo by phase
% idx_layer_n1                     = round(Idx_resin(end) - Idx_oneply / 2);
% idx_layer_n2                     = min(length(inph_ascan), round(Idx_resin(end) + Idx_oneply / 2));
% [pkt(end), lct(end)]             = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2) - (phase_0 - pi)));
% lct(end)                         = lct(end) + idx_layer_n1 - 1;

% fillna by average
% lct        = round(fx_inpaint_nans(lct, 1));

% % calculate errors
if length(lct)>length(Idx_resin)
    errors_pos = (lct(1:length(Idx_resin))-Idx_resin) / Idx_oneply;
else
    errors_pos = (lct-Idx_resin) / Idx_oneply;
end

% find the nearest one to calculate the errors
% errors_pos = NaN(1, length(Idx_resin));
% for idx_l = 1:length(Idx_resin)
%     [~, I_err]        = min(abs(lct - Idx_resin(idx_l)));
%     errors_pos(idx_l) = (lct(I_err) - Idx_resin(idx_l))/Idx_oneply;
% end

% % debug
% figure, plot(inph_ascan);
% hold on;
% scatter(idx_layer_n1s(idx_layer_n1s~=1), ...
%     inph_ascan(idx_layer_n1s(idx_layer_n1s~=1)), 'rd');
% hold on;
% scatter(idx_layer_n2s(idx_layer_n2s~=1), ...
%     inph_ascan(idx_layer_n2s(idx_layer_n2s~=1)), 'gs');
% hold on;
% scatter(lct(~isnan(lct)), inph_ascan(round(lct(~isnan(lct)))), 'yo', 'fill');

end

