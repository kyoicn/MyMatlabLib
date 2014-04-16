function re = intervaled(data, bin)
if (bin <= 0 || rem(bin, 1) ~= 0)
    error('invalid # of bins');
end

si = size(data);
re = zeros(bin, si(2));
if (si(1) * si(2) == 0) 
    return;
end

mi = min(data(:, 1));
ma = max(data(:, 1));
gap = (ma - mi) / bin;
[~, in] = sort(data(:, 1));
sd = data(in, :);

c = 0;
p = 1;
for i = 1 : bin
    while (p + c <= si(1) && sd(p + c) <= mi + i * gap) c = c + 1; end
    if (c > 0)
        re(i, :) = mean(sd(p:p+c-1, :), 1);
        p = p + c;
        c = 0;
    end
    if (p > si(1)) break; end
end