function [c_p, anglex, angley, anglez] = fx_decomSS_surface(ST, front_I, rear_I, type)
% form Structure Tensor in a 3D cell 
% calculate eigenvalues and eigenvectors

% choose the filter method by G_tau

if nargin==3
    type = 'planar';
else
    typestrings = {'planar', 'linear'};
    type = validatestring(type, typestrings);
end

switch type
    case 'planar'
        index_lambda = 1;
    case 'linear'
        index_lambda = 3;
end


[lx, ly, lz] = size(ST(:, :, :, 1, 1));

disp("size xyz");
disp([lx ly lz]);

% V      = cell(lx, ly, lz);
c_p    = single(nan(lx, ly, lz));
anglex = single(nan(lx, ly, lz));
angley = single(nan(lx, ly, lz));
anglez = single(nan(lx, ly, lz));

disp("start");

for i=1:lx
    for j=1:ly
        for k=  1:lz
            if isempty(front_I) || (k >= max(1, round(front_I(i, j))) && k <= min(lz, round(rear_I(i, j))))          
                matrix = squeeze(ST(i, j, k, :, :));
                [vector, d_d] = eig(matrix);
               % Extract the eigenvalues from the diagonal, then sort the resulting vector in ascending order. The second output from sort returns a permutation vector of indices.
                [d, ind]      = sort(diag(d_d)); 
    %             d = abs(d); % just extract the real part
    %             angles(i,j,k) = vector(ind(1), 3) / sqrt(vector(ind(1), 1)^2 + vector(ind(1), 2)^2);
                % disp(d_d);
                c_p(i, j, k) = (d(end) - d(end-1)) / d(end);
                vector1 = vector(:, ind(end)); % whose columns are the corresponding right eigenvectors,
                %calibrating the pixels scales
                %************ need more exploration here ! **********
%                 vector1(3) = vector1(3) / 0.012 * 0.02; % the 3rd indicate the z angle?
                anglex(i, j, k) = acos(vector1(1)) - pi/2;
                angley(i, j, k) = acos(vector1(2)) - pi/2;
                anglez(i, j, k) = acos(vector1(3)) - pi/2;
%                 V{i, j, k} = vector1;
            else
                c_p(i, j, k) = 0;
                anglex(i, j, k) = NaN;
                angley(i, j, k) = NaN;
                anglez(i, j, k) = NaN;
            end
        end
    end
    clc;
    disp(i + "/" + lx);
end

end


