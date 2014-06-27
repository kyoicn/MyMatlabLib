% Minimal coherence sampling
% samples = mc_sample(m, base)

function [samples, bestcohs, samples_sequence] = mc_sample_debug(m, b, rk, debug)
    % so far only deal with 2D
    if (nargin < 3) rk = 0; debug = 0;
    elseif (nargin < 4) debug = 0;
    end
    dim = size(b);
    n = dim(2);
    choose = max(1, floor(rk * n));
    if (n == 0) samples = []; bestcohs = 0; samples_sequence = [];
    elseif (n <= m) samples = ones(n, 1); bestcohs = coh(b); samples_sequence = samples;
    else
        samples = zeros(n, 1);
        samples_sequence = zeros(n, m);
        bestcohs = zeros(m, 1);
        for sam = 1 : m
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
            choose = randi(min(length(pos), floor(rk*n)));
            if (length(pos) >= choose) picked = pos(choose);
            else picked = pos(length(pos)); end
            bestcohs(sam) = cohs(picked);
            samples(idx(picked)) = 1;
            samples_sequence(:, sam) = samples;
            if (debug) fprintf('(%d,DCT4,%.1f):%d-th minial coh=%f, %d picked in round%d\n',...
                                m, rk, choose, bestcohs(sam), idx(picked), sam); end
        end
    end
end