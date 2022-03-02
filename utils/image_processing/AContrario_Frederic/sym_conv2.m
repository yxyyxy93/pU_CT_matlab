function C = sym_conv2(A,B)
% symmetrically padded convolution
% from http://www.mathworks.fr/support/solutions/en/data/1-308SBW/index.html?product=ML&solution=1-308SBW

xLayer=length(B);
A_x = padarray(A,[xLayer xLayer],'symmetric');
C_x = conv2(A_x,B,'same');
C = C_x(xLayer+1:end-xLayer, xLayer+1:end-xLayer);

end
