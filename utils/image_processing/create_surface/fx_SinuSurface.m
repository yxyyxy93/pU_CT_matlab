function surface = fx_SinuSurface(h, omega, phase1, phase2, x_mesh, y_mesh)
% create sinusoidal surface
% h: amplitdue
% x_mesh, y_mesh: x and y ranges

surface = h .* (sin(omega.*x_mesh + phase1) + sin(omega.*y_mesh + phase2));

end

