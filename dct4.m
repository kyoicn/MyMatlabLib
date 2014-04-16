%%DCT-IV

function [trans, basis] = dct4(origin)
dim = size(origin);
if (dim(1) >= dim(2))
    n = dim(1);
else
    n = dim(2);
    origin = origin';
end
basis = zeros(n);
for i = 0 : n - 1
    for j = 0 : n - 1
        basis(i + 1, j + 1) = cos(pi * (j + 0.5) * (i + 0.5) / n);
    end
end
basis = basis * sqrt(2 / n);
trans = basis * origin;