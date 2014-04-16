% 1-Dimension Discrete Cosine Tranform

function [trans, sparse, count, basis] = dct2(origin, ignore)
basisLength = length(origin);
count = 0;

for i = 0 : (basisLength - 1)
    for j = 0 : (basisLength - 1)
        basis(i + 1, j + 1) = cos(pi * (2 * j + 1) * i / 2 / basisLength);
    end
end

trans = basis * origin';
maxElt = abs(trans(1));
for i = 1 : basisLength
    if abs(trans(i)) > maxElt
        maxElt = abs(trans(i));
    end
end

for i = 1 : basisLength
    if abs(trans(i)) < ignore * maxElt
        sparse(i) = 0;
    else
        sparse(i) = trans(i);
        count = count + 1;
    end
end
sparse = sparse';