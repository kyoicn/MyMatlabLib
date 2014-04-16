%%DCT-IV

function basis = base_dct4ex(m, n)
if (m < n) d1 = n; % d1 is the larger dimension
else d1 = m; end
d2 = m + n - d1;

basis = zeros(d1, d2);
for i = 0 : d1 - 1
    for j = 0 : d2 - 1
        basis(i + 1, j + 1) = cos(pi * (j + 0.5) * (i + 0.5) / d2);
    end
end
basis = basis * sqrt(2 / n);
if (m < n) basis = basis'; end