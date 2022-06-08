function h = fx_handwritten_3DGabor(thetas, wl, SFB, SAR)
       % sigma
    sigma  = wl/pi*sqrt(log(2)/2)*(2^SFB+1)/(2^SFB-1);
    SigmaX = sigma;
    SigmaY = sigma ./ SAR;
    SigmaZ = sigma * SAR;
    
    % SpatialKernel
    rx = ceil(7*SigmaX);
    ry = ceil(7*SigmaY);
    rz = ceil(7*SigmaZ);
    
    r     = max([rx, ry, rz]);
    [X, Y, Z] = meshgrid(-r:r, -r:r, -r:r);
%     kSize = [2*r+1, 2*r+1, 2*r+1];

    % Gaussain modulation
    phi = 0;
%     X_prime = X .* cosd(orientation) - Y .*sind(orientation);
%     Y_prime = X .* sind(orientation) + Y .*cosd(orientation);
%     Z_prime = Z;
	R = rotation(thetas);
    Z_prime = Z * R(1,1) + Y * R(1,2) + X * R(1,3);
    Y_prime = Z * R(2,1) + Y * R(2,2) + X * R(2,3);
    X_prime = Z * R(3,1) + Y * R(3,2) + X * R(3,3);
    
    hGaussian  = exp( -1/2*( X_prime.^2 ./ SigmaX^2 + Y_prime.^2 ./ SigmaY^2 + Z_prime.^2 ./ SigmaZ^2));
%     hGaussian  = exp( -1/2*( X_prime.^2 ./ SigmaX^2 + Y_prime.^2 ./ SigmaY^2));
    hGaborEven = hGaussian.*cos(2*pi.*X_prime ./ wl+phi);
    hGaborOdd  = hGaussian.*sin(2*pi.*X_prime ./ wl+phi);

    %
%     hGaborEven(Z_prime > 0) = 0;
%     hGaborOdd(Z_prime < 0)  = 0;
%     
    h = complex(hGaborEven,hGaborOdd);
end

function R = rotation(theta)
    R_x = [1,         0,             0;
        0,         cosd(theta(1)), -sind(theta(1));
        0,         sind(theta(1)), cosd(theta(1))];
         
    R_y = [cosd(theta(2)),  0,      sind(theta(2));
        0,               1,      0;
        -sind(theta(2)), 0,      cosd(theta(2))];
                 
    R_z = [cosd(theta(3)),    -sind(theta(3)),    0;
        sind(theta(3)),    cosd(theta(3)),     0;
        0,                 0,                  1];
        
    R = R_z * R_x * R_y;
end