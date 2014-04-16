%%DCT-IV

function basis = base_dct4(n)
basis = zeros(n);
for i = 0 : n - 1
    for j = 0 : n - 1
        basis(i + 1, j + 1) = cos(pi * (j + 0.5) * (i + 0.5) / n);
    end
end
basis = basis * sqrt(2 / n);