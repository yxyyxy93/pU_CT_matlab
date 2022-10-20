function u_n = fx_mumfordshah(f, alpha, lambda, eps)
% Summary of this function goes here
% Detailed explanation goes here
% f: iamge 
% alpha, lambda, eps: algorithm para

[m, n, k] = size(f); % k = channels, 1 for gray scalar

u_n = f;
umac_n = u_n;
p_nx = 0;
p_ny = 0;
tau = 1 / (2*2);
sigma = 1 / 2;


%% solver
while true
    % dual ascent in p
    fx = umac_n(:, 2:end, :) - umac_n(:, 1:end-1, :);
    fx = padarray(fx, [0 1 0], 'post');
    fy = umac_n(2:end, :, :) - umac_n(1:end-1, :, :);
    fy = padarray(fy, [1 0 0], 'post');
    
    p_tilx = p_nx + sigma * fx;
    p_tily = p_ny + sigma * fy;
    
    p_nx = 2*alpha / (sigma + 2*alpha) .* p_tilx; % update to p^(n+1)
    p_ny = 2*alpha / (sigma + 2*alpha) .* p_tily; % update to p^(n+1)
    
    p_nx(abs(p_tilx) > sqrt(lambda / alpha * sigma * (sigma + 2*alpha))) = 0;
    p_ny(abs(p_tily) > sqrt(lambda / alpha * sigma * (sigma + 2*alpha))) = 0;

    % primal descent in u

    div = padarray(p_nx(:, 2:end, :) - p_nx(:, 1:end-1, :), [0 1 0], 'pre') ...
        + padarray(p_ny(2:end, :, :) - p_ny(1:end-1, :, :), [1 0 0], 'pre');
    
    u_til = u_n + tau * div;
    u_pre = u_n;
    u_n =  (u_til + 2*tau*f) / (1 + 2*tau); % update to u^(n+1)
    
    % extrapolation step
    theta = 1 / sqrt(1 + 4*tau);
    tau = theta * tau;
    sigma = sigma / theta;
    umac_n = u_n + theta * (u_n - u_pre);
    
    % exit condition 
    increase = sum((u_n - u_pre).^2, 'all')^(1/2);
    disp(increase);
    if increase  < eps
        break;
    end
end


end
