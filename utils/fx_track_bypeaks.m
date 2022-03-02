function [lct, pkt, errors_pos, inam_ascan] = fx_track_bypeaks(Idx_resin, ori_signal, threshold, power_noise)
% find the peaks from the original signal by a rule
% get the ply thickness

P                 = polyfit(1:length(Idx_resin), Idx_resin, 1);
Idx_oneply        = P(1);

% find the peaks in inam, as the front-wall
if isempty(find(ori_signal<0, 1)) % inst. amp envelop
    inam_ascan        = ori_signal;
else
    inam_ascan    = abs(hilbert(ori_signal));
end
MinPeakHeight     = threshold * max(inam_ascan);
% MinPeakDistance        = length(inam_ascan) / 2;
[~, front_I]      = max(inam_ascan);
C_max             = max(ori_signal);

% back_I                 = round(Idx_walls(end));
%do a parabolic estimate
[lct(1), pkt(1), ~]    = crit_interp_p(inam_ascan(front_I+(-1:1)), front_I+(-1:1));
lct(1)                 = round(lct(1)); % get the int

% [back_I, ~, ~]         = crit_interp_p(inam_ascan(back_I+(-1:1)), back_I+(-1:1));

pre_lct           = lct(1);

for n_layer = 2:length(Idx_resin) % from 1 for simulation, from 2 for exp. 
    % idx_layer_n1                         = round(Idx_resin(n_layer)-Idx_oneply/2);
    % idx_layer_n2                         = round(Idx_resin(n_layer)+Idx_oneply/2);
    idx_layer_n1 = min(length(inam_ascan)-3, round(pre_lct + Idx_oneply/2));
    idx_layer_n2 = min(length(inam_ascan), round(pre_lct + Idx_oneply*3/2));
    idx_layer_n2 = max(idx_layer_n1+3, idx_layer_n2); % Data set must contain at least 3 samples.
    inam_oneply  = inam_ascan(idx_layer_n1: idx_layer_n2);
    %     ascan_oneply = ascan(idx_layer_n1: idx_layer_n2);
    %     adiff_oneply = diff(inam_oneply)/mean(inam_oneply);
    %         find local maximum
    [pks_n, locs_n]  = findpeaks(inam_oneply, 'MinPeakHeight', MinPeakHeight);
    %     %     peakcoe_sig      = (max(pks_n)-mean(inam_oneply))^2 / mean(inam_oneply);
    %     peakcoe_sig      = (max(pks_n) - mean(inam_oneply) ) / C_max;
    %     %     peakcoe_sig      = max(pks_n) / mean(inam_oneply);
    %     if peakcoe_sig < 1e-10 * power_noise % the kurtosis is less than that from noise
    %         %         break; % just stop
    %         lct(n_layer) = NaN;
    %         pkt(n_layer) = NaN;
    %         pre_lct      = pre_lct + Idx_oneply; % updata pre_lct
    %         %         disp(['nan']);
    %     else
    if isempty(locs_n)
        %         % detrend and find new maximum
        %         inam_oneply_detrend       = detrend(inam_oneply);
        pos                       = randi(length(inam_oneply));
        lct_local                 = pos;
        lct(n_layer)              = lct_local + idx_layer_n1 - 1;
        pre_lct                   = lct(n_layer); % updata pre_lct
    else
        [pkt(n_layer), lct_local] = max(pks_n);
        lct(n_layer)              = locs_n(lct_local) + idx_layer_n1 - 1;
        pre_lct                   = lct(n_layer); % updata pre_lct
    end
    %     end
    %     % check the replicated track points.
    %     if lct(n_layer)==lct(n_layer-1)
    %         lct(n_layer)             = NaN;
    %     elseif lct(n_layer)>= Idx_resin(end) % exceed the rear surface
    %         lct(n_layer)             = NaN;
    %         break;
    %     end
    % %     by maximum
    %     if (idx_layer_n2 <= length(inam))
    %         [pkt(n_layer), lct_local]     = max(inam(idx_layer_n1: idx_layer_n2));
    %         lct(n_layer)                  = lct_local + idx_layer_n1 - 1;
    %     else
    %         break;
    %     end
end
% calculate errors
errors_pos = (lct-Idx_resin)/Idx_oneply;

end

