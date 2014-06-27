% Min-Max distance Pair-wise sampling
% samples = pw_sample(m, cost)

function samples = pw_sample(m, cost)
    dim = size(cost);
    n = prod(dim);
    
    if (n == 0) samples = [];
    elseif (n < m) samples = ones(dim);
    else
        samples = zeros(dim);
        m0 = m;
        if (2*m0 > n) m = n - m; end
        cands = randperm(n, 2*m);
        % wait to generalize to higher dimensions
        coy = ceil(cands./dim(1));
        cox = cands - (coy-1)*dim(1);
        [tx1 tx2] = meshgrid(cox, cox);
        [ty1 ty2] = meshgrid(coy, coy);
        dis = (tx2-tx1).^2 + (ty2-ty1).^2 + diag(inf*ones(1, 2*m));
        for sam = 1 : m
            [mins,mxs] = min(dis);
            [~,myi] = min(mins);
            mxi = mxs(myi);
            p1 = [cox(mxi) coy(mxi)];
            p2 = [cox(myi) coy(myi)];
            if (cost(p1(1), p1(2)) < cost(p2(1), p2(2)))
                if (2*m0 > n) samples(p2(1), p2(2)) = 1;
                else samples(p1(1), p1(2)) = 1; end
            else
                if (2*m0 > n) samples(p1(1), p1(2)) = 1;
                else samples(p2(1), p2(2)) = 1; end
            end
            dis([mxi myi], :) = inf*ones(2, 2*m);
            dis(:, [mxi myi]) = inf*ones(2*m, 2);
        end
        if (2*m0 > n) samples = ones(dim)-samples; end
    end
end