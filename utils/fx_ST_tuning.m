function [c_p, anglex, angley, anglez] = fx_ST_tuning(inph, grad_kernel, ...
    s_kernel, integ_scale, front_I, rear_I, xidx, yidx)
%%
tic;
% time record
timeElapsed = toc;
disp(['deconve: ', num2str(timeElapsed)]);

% 3D Gaussian filter
inph_cos = cos(inph);
inph_sin = sin(inph);

% smoothing scale
sigma1 = s_kernel;

inph_cos_filter = imgaussfilt3(inph_cos, sigma1, 'padding', 'replicate');
inph_sin_filter = imgaussfilt3(inph_sin, sigma1, 'padding', 'replicate');

% fx_showBS(2, inph_ex, num_bs, step_bs, front_I, rear_I_PhaseDrived, B_type);
% fx_showBS(3, inph_ex_filter, num_bs, step_bs, front_I, rear_I_PhaseDrived, B_type);

% calculation of gradients
% inph_ex_filter = inph_ex; % no filter
[Gx_cos, Gy_cos, Gz_cos] = fx_imgradientxyz(inph_cos_filter, grad_kernel);
[Gx_sin, Gy_sin, Gz_sin] = fx_imgradientxyz(inph_sin_filter, grad_kernel);

Gx = inph_cos.*Gx_sin - inph_sin.*Gx_cos;
Gy = inph_cos.*Gy_sin - inph_sin.*Gy_cos;
Gz = inph_cos.*Gz_sin - inph_sin.*Gz_cos;


% form the structure-tensor
ST = fx_formStructureTensor(Gx, Gy, Gz);
ST_in = zeros(size(ST));

% integration scale
sigma2 = integ_scale;
for i=1:3
    for j=1:3
        ST_in(:, :, :, i, j) = imgaussfilt3(ST(:, :, :, i, j), sigma2);
    end
end

% time record
timeElapsed = toc;
disp(['form_tensor: ', num2str(timeElapsed)]);

% x section
ST_in_x = ST_in(xidx, :, :, :, :);
front_I_temp = front_I(xidx, :);
rear_I_temp = rear_I(xidx, :);
[~, c_p, anglex, ~, ~] = fx_decomSS_surface(ST_in_x, front_I_temp, rear_I_temp);

% x section
ST_in_y = ST_in(:, yidx, :, :, :);
front_I_temp = front_I(:, yidx);
rear_I_temp = rear_I(:, yidx);
[~, c_p, ~, angley, anglez] = fx_decomSS_surface(ST_in_y, front_I_temp, rear_I_temp);

% time record
timeElapsed = toc;
disp(['extraxt angles: ', num2str(timeElapsed)]);


end

