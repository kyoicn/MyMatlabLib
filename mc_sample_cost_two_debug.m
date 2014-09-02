% Minimal coherence sampling
% samples = mc_sample(m, base)

function [samples, bestcohs, samples_sequence] = mc_sample_cost_two_debug(m, b, cost, seed, debug)
    % so far only deal with 2D
    if (nargin < 4) seed = 0; debug = 0;
    elseif (nargin < 5) debug = 0;
    end
    dim = size(b);
    n = dim(2);
    sumcost = sum(cost);
    acost = 0;
    if (n == 0) samples = []; bestcohs = 0; samples_sequence = [];
    elseif (n <= m) samples = ones(n, 1); bestcohs = coh(b); samples_sequence = samples;
    else
        samples = zeros(n, 1);
        samples_sequence = zeros(n, m);
        bestcohs = zeros(m, 1);
        for sam = 1 : m
            if (sam == 1)
                if (seed == 0) [~, seed] = min(cost); end
                c = cost(seed);
                samples(seed) = 1;
                samples_sequence(:,sam) = samples;
                bestcohs(sam) = 1;
                acost = acost + c;
                if (debug) fprintf('(%d,min-coh-two):seed %d, coh=1, cost=%f%%\n', m, seed, acost/sumcost*100); end
                continue;
            end
            cohs = zeros(1, n-sam+1);
            samplex = find(samples > 0);
            idx = find(samples == 0);
            for candi = 1 : length(idx)
                cand = idx(candi);
                part = b([samplex; cand], :);
                cohs(candi) = coh(part);
            end
%             [~, pos] = min(cohs);
            [~, pos] = sort(cohs);
            [c, choose] = min(cost(idx(pos(1:min(length(pos), 2)))));
            picked = pos(choose);
            acost = acost + c;
            bestcohs(sam) = cohs(picked);
            samples(idx(picked)) = 1;
            samples_sequence(:, sam) = samples;
            if (debug) fprintf('(%d,min-coh-two):%d-th minial coh=%f, cost=%f%%, %d picked in round%d\n',...
                                m, choose, bestcohs(sam), acost/sumcost*100, idx(picked), sam); end
        end
    end
end