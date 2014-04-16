function [ b, V ] = pca( A )
% [b,v] = pca(A)
% A = V*b;
CorrA = A*A';
[S,V,~] = svd(CorrA);
b = S'*A;

end

