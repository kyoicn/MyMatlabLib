% Minimal coherence sampling
% samples = mc_sample(m, base)

function [samples, bestcohs, samples_sequence] = mc_sample_cost_debug(m, b, cost, rk, seed, debug)
    % so far only deal with 2D
    if (nargin < 4) rk = 0; seed = 0; debug = 0;
    elseif (nargin < 5) seed = 0; debug = 0;
    elseif (nargin < 6) debug = 0;
    end
    dim = size(b);
    n = dim(2);
    sumcost = sum(cost);
    acost = 0;
    if (n == 0) samples = []; bestcohs = 0; samples_sequence = [];
    elseif (n < m) samples = 1:n; bestcohs = coh(b); samples_sequence = zeros(n,1);
    else
        samples = zeros(m, 1);
        flag = zeros(n, 1);
        samples_sequence = zeros(n, m);
        bestcohs = zeros(m, 1);
        for sam = 1 : m
            if (sam == 1)
                if (seed == 0) [~, seed] = min(cost); end
                c = cost(seed);
                samples(sam) = seed;
                flag(seed) = 1;
                samples_sequence(:,sam) = flag;
                bestcohs(sam) = 1;
                acost = acost + c;
                if (debug) fprintf('(%d,%d-min):seed %d, coh=1, cost=%f%%\n', m, floor(rk*n), seed, acost/sumcost*100); end
                continue;
            end
            cohs = zeros(1, n-sam+1);
            samplex = find(flag > 0);
            idx = find(flag == 0);
            for candi = 1 : length(idx)
                cand = idx(candi);
                part = b([samplex; cand], :);
                cohs(candi) = coh(part);
            end
%             [~, pos] = min(cohs);
            [~, pos] = sort(cohs,'descend');
            [c, choose] = min(cost(idx(pos(1:min(length(pos), floor(rk*n))))));
            picked = pos(choose);
            acost = acost + c;
            bestcohs(sam) = cohs(picked);
            samples(sam) = idx(picked);
            flag(idx(picked)) = 1;
            samples_sequence(:, sam) = flag;
            if (debug) fprintf('(%d,%d-min):%d-th minial coh=%f, cost=%f%%, %d picked in round%d\n',...
                                m, floor(rk*n), choose, bestcohs(sam), acost/sumcost*100, idx(picked), sam); end
        end
    end
end