% Minimal norm sampling: with cost
% samples = minnorm_sample(m, base)
%
% m: number of samples
% b: base
% id: identity column index
% rk:

function [samples, bestnorm, samples_sequence] = minnorm_sample_cost(m, b, id, cost, rk, seed, random, debug)
    % so far only deal with 2D
    if (nargin < 5) rk = 0; seed = 0; random = 0; debug = 0;
    elseif (nargin < 6) seed = 0; random = 0; debug = 0;
    elseif (nargin < 7) random = 0; debug = 0;
    elseif (nargin < 8) debug = 0;
    end
    dim = size(b);
    n = dim(2);
    sumcost = sum(cost);
    acost = 0;
    cols = 1:n;
    cols(id) = [];
    if (n == 0) samples = []; bestnorm = 0; samples_sequence = [];
    elseif (n < m) samples = 1:n; bestnorm = mcs(b(:,cols)); samples_sequence = ones(n,1);
    else
        samples = zeros(m, 1);
        flag = zeros(n, 1);
        samples_sequence = zeros(n, m);
        bestnorm = zeros(m, 1);
        for sam = 1 : m
            if (sam == 1)
                if (seed == 0 && random ~= 0) seed = randi(n);
                elseif (seed == 0 && random == 0) [~, seed] = min(cost); end
                samples(sam) = seed;
                c = cost(seed);
                acost = acost + c;
                flag(seed) = 1;
                samples_sequence(:,sam) = flag;
                bestnorm(sam) = 1;
                if (debug) fprintf('(%d,B,%.1f):seed %d, norm=1, cost=%.1f%%(%.1f%%)\n',...
                        m, rk, seed, acost/sumcost*100, (acost/sumcost-sam/n)*100); end
                continue;
            end
            norms = zeros(1, n-sam+1);
            samplex = find(flag > 0);
            idx = find(flag == 0);
            for candi = 1 : length(idx)
                cand = idx(candi);
                part = b([samplex; cand], cols);
                norms(candi) = mcs(part);
            end
%             [~, pos] = min(cohs);
            [~, pos] = sort(norms);
            [c, choose] = min(cost(idx(pos(1:min(length(pos), floor(rk*n))))));
            %if (random ~= 0) choose = randi(min(length(pos), floor(rk*n)));end
            %if (length(pos) < choose) choose = length(pos); end
            picked = pos(choose);
            bestnorm(sam) = norms(picked);
            samples(sam) = idx(picked);
            acost = acost + c;
            flag(idx(picked)) = 1;
            samples_sequence(:, sam) = flag;
            if (debug) fprintf('(%d,B,%.1f):%d-th minial norm=%f, cost=%.1f%%(%.1f%%), %d picked in round%d\n',...
                                m, rk, choose, bestnorm(sam), acost/sumcost*100, (acost/sumcost-sam/n)*100, idx(picked), sam); end
        end
    end
end