function re = linearspline(n, sample, index)
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
if (in(1) > 1)
    xl = re(in(1)) - re(in(1) + 1);
    for dc = in(1) - 1 : -1 : 1
        re(dc) = re(dc + 1) + xl;
    end
end
if (in(m) < n)
    xl = re(in(m)) - re(in(m) - 1);
    for dc = in(m) + 1 : 1 : n
        re(dc) = re(dc - 1) + xl;
    end
end