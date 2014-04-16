%% q-order chebyshev polynomial base

function [trans, base] = cpb(origin, q)
dim = size(origin);
if (dim(1) >= dim(2))
    n = dim(1);
else
    n = dim(2);
    origin = origin';
end
base = zeros(n);

for i = 1 : n
    for j = 1 : q
        base(i, j) = cos(j * acos((2*i) / (n+1) - 1));
        base(i, j) = base(i, j) / sqrt(1 - ((2*i) / (n+1) - 1)^2);
    end
end

trans = pinv(base) * origin;