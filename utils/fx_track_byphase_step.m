function [lct, pkt, errors_pos, inam_ascan, inph_ascan] = fx_track_byphase_step(Idx_resin, ori_signal, fs, threshold, power_noise)

% calculate average TOF ply
P          = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply = P(1);
%
inph_ascan = angle(hilbert(ori_signal));
inam_ascan = abs(hilbert(ori_signal));
infq_ascan = 1/(2*pi) * fs * diff(unwrap(inph_ascan));

%
front_I = round(Idx_resin(1));
back_I  = round(Idx_resin(end));
lct      = NaN(1, length(Idx_resin));
lct(1)   = front_I;
lct(end) = back_I;
phase_0 = interp1(round(front_I)+(-1:1), inph_ascan(round(front_I)+(-1:1)), front_I);

% find the peaks in inam, as the front- and back- walls
% MinPeakHeight          = threshold * max(inam_ascan);
% MinPeakDistance        = length(inam_ascan) / 2;
% [pks_walls, Idx_walls] = findpeaks(inam_ascan, 1:length(inam_ascan), 'MinPeakHeight',  MinPeakHeight, 'MinPeakDistance', MinPeakDistance);
% if isempty(pks_walls)
%     front_I = nan;
%     lct     = [];
%     pkt     = [];
%     errors_pos = [];
%     Idx_resin = []; % just end the function.
% else
%     front_I = round(Idx_walls(1));
%     back_I  = round(Idx_walls(end));
%     % set the region of each layer and find the postion where
%     % the phase is closest to - pi / 2
%     pkt                     = NaN(1, length(Idx_resin));
%     lct                     = NaN(1, length(Idx_resin));
%     lct(1)                  = front_I;
%     lct(end)                = back_I;
%     pkt(1)                  = pks_walls(1);
%     pkt(end)                = pks_walls(end);
%     % interpolation
%     phase_0                 = interp1(round(front_I)+(-1:1), inph_ascan(round(front_I)+(-1:1)), front_I);
%     pkt(1)                  = phase_0; % all the pkt are phases!
% end

% % find the back-wall echo by phase at first
% idx_layer_n1            = max(front_I, round(Idx_resin(end) - Idx_oneply / 2));
% idx_layer_n2            = min(length(inph_ascan), round(Idx_resin(end) + Idx_oneply / 2));
% if phase_0 > 0
%     phase_back = phase_0 - pi;
% else
%     phase_back = phase_0 - pi + 2*pi;
% end
% [pkt(end), lct_back] = min(abs(inph_ascan(idx_layer_n1: idx_layer_n2)-phase_back));
% lct(end)                = lct_back + idx_layer_n1 - 1;
% pkt(end)                = inph_ascan(round(lct(end))); % all the pkt are phases!
pre_lct                 = front_I;

for n_layer = 2:length(Idx_resin)-1 % from 1 for simulation, from 2 for exp. until last interply before rear surface
    % idx_layer_n1                         = round(Idx_resin(n_layer)-Idx_oneply/2);
    % idx_layer_n2                         = round(Idx_resin(n_layer)+Idx_oneply/2);
    idx_layer_n1      = min(length(infq_ascan)-2, round(pre_lct + Idx_oneply/2));
    idx_layer_n2      = min(length(infq_ascan), round(pre_lct + Idx_oneply*3/2));
    if idx_layer_n2 >= lct(end) % exceed the back surface bound
        break;
    end
    % ascan_oneply = ascan(idx_layer_n1: idx_layer_n2);
    % adiff_oneply = diff(inam_oneply)/mean(inam_oneply);
    % find local maximum
    % peakcoe_sig      = (max(pks_n)-mean(inam_oneply))^2 / mean(inam_oneply);
    infq_ascan_oneply = infq_ascan(idx_layer_n1: idx_layer_n2);
    peakcoe_sig       = std(infq_ascan_oneply);
    if peakcoe_sig > power_noise/2  * 1e10
        lct(n_layer) = NaN;
        pkt(n_layer) = NaN;
        pre_lct      = pre_lct + Idx_oneply; % updata pre_lct
    else
        % find one nearest point and apply interpolation
        inph_ascan_oneply         = inph_ascan(idx_layer_n1: idx_layer_n2);
        %in order to only find positive numbers
        inph_check                = inph_ascan_oneply - (phase_0-pi/2);
        if sum(inph_check>0)>0 % take the case of all "NaN" into account
            inph_check(inph_check<=0) = nan; %replace negative numbers and the zero with nan
            [~, I_minpos]             = min(inph_check, [], 'omitnan'); %find values
        else
            I_minpos = NaN;
        end
        %in order to only find negative numbers
        inph_check                = inph_ascan_oneply - (phase_0-pi/2);
        if sum(inph_check<=0)>0 % take the case of all "NaN" into account
            inph_check(inph_check>0) = nan; %replace positive numbers and the zero with nan
            [~, I_maxneg]            = max(inph_check, [], 'omitnan'); %find values
        else
            I_maxneg = NaN;
        end
        nearest_ps                    = [I_maxneg I_minpos];
        nearest_ps(isnan(nearest_ps)) = [];               
        %     % check whether it exceeds the bound
        %     left_b               = max(1, nearest_p-1);
        %     right_b              = min(length(inph_ascan_oneply), nearest_p+1);
        %     nearest_ps           = left_b:1:right_b;
        x                    = inph_ascan_oneply(nearest_ps);
        [~, ind]             = unique(x); % ind = index of first occurrence of a repeated value
        if length(ind)==1
            lct(n_layer) = nearest_ps(ind)+idx_layer_n1-1;
            pkt(n_layer) = phase_0-pi/2;
            pre_lct      = lct(n_layer); % updata pre_lct
        elseif isempty(ind)
            lct(n_layer) = NaN;
            pkt(n_layer) = NaN;
            pre_lct      = pre_lct + Idx_oneply; % updata pre_lct
        else
            vq1          = interp1(x(ind), nearest_ps(ind), phase_0-pi/2, 'linear');
            lct(n_layer) = vq1+idx_layer_n1-1;
            pkt(n_layer) = phase_0-pi/2;
            pre_lct      = vq1+idx_layer_n1-1; % updata pre_lct
        end
    end
end

% calculate errors
errors_pos = (lct - Idx_resin)/Idx_oneply;

end

