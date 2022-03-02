function [c_p, c_s, c_l] = fx_decomSS_checkfeature(ST)
% form Structure Tensor in a 3D cell 
% calculate eigenvalues and eigenvectors

% choose the filter method by G_tau

[lx, ly, lz] = size(ST(:, :, :, 1, 1));

c_p = cell(lx, ly, lz);
c_s = cell(lx, ly, lz);
c_l = cell(lx, ly, lz);

disp("start");

for i=1:lx
    for j=1:ly
        for k=1:lz

            matrix = squeeze(ST(i,j,k,:,:));
            [~, d_d] = eig(matrix);
            
           % Extract the eigenvalues from the diagonal, then sort the resulting vector in ascending order. 
           % The second output from sort returns a permutation vector of indices.
            [d, ~]  = sort(diag(d_d)); 
            
            %check the feature of tensor 
            c_p{i,j,k} = (d(end)-d(end-1))/d(end);
            c_s{i,j,k} = d(1)/d(end);
            c_l{i,j,k} = (d(2)-d(1))/d(end);
       
        end
    end
    disp(i + "/" + lx);
end

end


