function re = linere(n, sample, index)
re = zeros(1, n);
si = size(sample);
if (si(1) == 1)
    sa = sample;
elseif (si(2) == 1)
    sa = sample';
end
si = size(index);
if (si(1) == 1)
    in = index;
elseif (si(2) == 1)
    in = index';
end
[in, tin] = sort(in, 2, 'ascend');
sa = sa(tin);
m = length(sa);

for j = 1 : m
	re(in) = sa;
	if (j > 1)
        gap = (sa(j) - sa(j - 1)) / (in(j) - in(j - 1));
        for q = in(j - 1) + 1 : in(j)
            re(q) = re(q - 1) + gap;
        end
	end
end