function A = fx_ThreeDimShift(A, offsets)
%Like circshift, but shifts are not circulant. Missing data are filled with zeros.
% shift last dimention in 3d dataset A 
%  B=noncircshift(A,offsets)

siz=size(A);

if size(offsets)~=siz(1:2)
   error("wrong size of offsets");
end

% B = zeros(siz);
for i=1:siz(1)
      for j=1:siz(2)
          A(i,j,:) = circshift(A(i,j,:), -round(offsets(i,j)));
      end
end

end