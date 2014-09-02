%% Local least cost sample

function samples = llc_sample(m, cost, r)
    dim = size(cost);
%     ndim = length(dim);
    n = prod(dim);
    if (n == 0) samples = []; return;
    elseif (n <= m) samples = ones(dim); return;
    end
    
    if (nargin < 3) r = max(1, sqrt(n/m)/2); end % adaptive radius
    samples = zeros(dim);
    bak_cost = cost;
    seeds = sort(randperm(n, m));
    % so far only for 2d
    coy = ceil(seeds./dim(1));
    cox = seeds - (coy-1)*dim(1);
    for seedi = 1 : m
        sx = cox(seedi);
        sy = coy(seedi);
        startx = max(1, floor(sx - r));
        starty = max(1, floor(sy - r));
        endx = min(dim(1), ceil(sx + r));
        endy = min(dim(2), ceil(sy + r));
        pcost = bak_cost(startx:endx, starty:endy);
        [tx ty] = meshgrid(1:endx-startx+1,1:endy-starty+1);
        tx = tx'; ty = ty';
        basep = [sx+1-startx sy+1-starty];
        basex = basep(1)*ones(size(pcost));
        basey = basep(2)*ones(size(pcost));
        dis = (tx-basex).^2+(ty-basey).^2;
        idx = find(dis(:) <= r^2);
        [~, offset_i] = min(pcost(idx));
        offset = idx(offset_i);
        yy = ceil(offset/size(pcost, 1));
        xx = offset - (yy-1) * size(pcost, 1);
        xx = startx + xx - 1;
        yy = starty + yy - 1;
        samples(xx,yy) = 1;
        bak_cost(xx, yy) = inf;
    end
end