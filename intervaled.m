function re = intervaled(data, bin, sortidx)
% group data into equally intervaled bins
% re = intervaled(data, bin, sortidx)
% data: input matrix
% bin: # of bins (>0)
% sortidx: index of reference column for sorting (default=1)

if (bin <= 0 || rem(bin, 1) ~= 0)
    error('invalid # of bins');
end
if (nargin < 3) sortidx = 1; end

si = size(data);

if (si(1) * si(2) == 0) 
    error('empty data');
%     return;
end
if (sortidx > si(2))
    error('sorting index exceeds index boundaries');
end

mi = min(data(:, sortidx));
ma = max(data(:, sortidx));
if (isinf(mi) || isinf(ma))
    error('infinite element exists in sorting column');
end
if (mi == ma && bin > 1)
    warning('min/max values in sorting column are identical, and the # of bins is greater than 1');
end

gap = (ma - mi) / bin;

re = zeros(bin, si(2));
counter = zeros(1, bin);
for i = 1 : si(1)
    cbin = floor((data(i, sortidx) - mi)/gap) + 1;
    if (cbin == bin+1) cbin = bin; end
    re(cbin, :) = re(cbin, :) + data(i, :);
    counter(cbin) = counter(cbin) + 1;
end
for i = 1 : bin
    if (counter(i) == 0) re(i, :) = nan;
    else re(i, :) = re(i, :) / counter(i); end
end