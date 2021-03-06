% Minimal coherence sampling
% samples = mc_sample(m, base)

function [samples, bestcohs, samples_sequence] = mc_sample_debug(m, b, rk, seed, random, debug)
    % so far only deal with 2D
    if (nargin < 3) rk = 0; seed = 0; random = 0; debug = 0;
    elseif (nargin < 4) seed = 0; random = 0; debug = 0;
    elseif (nargin < 5) random = 0; debug = 0;
    elseif (nargin < 6) debug = 0;
    end
    dim = size(b);
    n = dim(2);
    if (random == 0) choose = max(1, floor(rk * n)); end
    if (n == 0) samples = []; bestcohs = 0; samples_sequence = [];
    elseif (n < m) samples = 1:n; bestcohs = coh(b); samples_sequence = ones(n,1);
    else
        samples = zeros(m, 1);
        flag = zeros(n, 1);
        samples_sequence = zeros(n, m);
        bestcohs = zeros(m, 1);
        for sam = 1 : m
            if (sam == 1)
                if (seed == 0) seed = randi(n); end
                samples(sam) = seed;
                flag(seed) = 1;
                samples_sequence(:,sam) = flag;
                bestcohs(sam) = 1;
                if (debug) fprintf('(%d,DCT4,%.1f):seed %d, coh=1\n', m, rk, seed); end
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
            [~, pos] = sort(cohs);
            if (random ~= 0) choose = randi(min(length(pos), floor(rk*n)));end
            if (length(pos) < choose) choose = length(pos); end
            picked = pos(choose);
            bestcohs(sam) = cohs(picked);
            samples(sam) = idx(picked);
            flag(idx(picked)) = 1;
            samples_sequence(:, sam) = flag;
            if (debug) fprintf('(%d,DCT4,%.1f):%d-th minial coh=%f, %d picked in round%d\n',...
                                m, rk, choose, bestcohs(sam), idx(picked), sam); end
        end
    end
end