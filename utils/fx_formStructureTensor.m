function ST = fx_formStructureTensor(Gx, Gy, Gz)
% form Structure Tensor in a 3D cell 

[lx, ly, lz] = size(Gx);

ST = zeros(lx, ly, lz, 3, 3);

disp("start forming the ST");

ST(:, :, :, 1, 1) = Gx.*Gx;
ST(:, :, :, 1, 2) = Gx.*Gy;
ST(:, :, :, 1, 3) = Gx.*Gz;

ST(:, :, :, 2, 1) = Gy.*Gx;
ST(:, :, :, 2, 2) = Gy.*Gy;
ST(:, :, :, 2, 3) = Gy.*Gz;

ST(:, :, :, 3, 1) = Gz.*Gx;
ST(:, :, :, 3, 2) = Gz.*Gy;
ST(:, :, :, 3, 3) = Gz.*Gz;

disp("ST formed");

end


