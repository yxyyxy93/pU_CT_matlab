function outmatrix=fx_3Dsubmatrix2matrix(inputCell, varargin)
%MAT2TILES - merge a cell array of sub-matrix to a big marix sizes
%

% 24/12/2019
% Xiaoyu Yang, PhD candidate in ZJU & UGent

[lcx, lcy] = size(inputCell);

outmatrix_row = inputCell{1,1};

for i=1:lcx   
    outmatrix_row = inputCell{i,1};
    for j=2:lcy
        outmatrix_row = [outmatrix_row inputCell{i,j}];
    end
    if (i==1)
        outmatrix = outmatrix_row;
    else
        outmatrix = [outmatrix; outmatrix_row];
    end
    disp(i);
end%loop over the dimensions
    
    
end
