% Minimal norm sampling
% samples = minnorm_sample(m, base)
%
% m: number of samples
% b: base
% id: identity column index
% rk: <= rk-th smallest norm value are considered 
% seed: initially selected row
% random: Obsoleted
% debug: display debug info

function [samples, bestnorm, samples_sequence] = minnorm_sample_test(m, b, id, rk, seed, random, debug)
    % so far only deal with 2D
    if (nargin < 4) rk = 0; seed = 0; random = 0; debug = 0;
    elseif (nargin < 5) seed = 0; random = 0; debug = 0;
    elseif (nargin < 6) random = 0; debug = 0;
    elseif (nargin < 7) debug = 0;
    end
    dim = size(b);
    n = dim(2);
    cols = 1:n;
    cols(id) = [];
%     if (random == 0) choose = max(1, floor(rk * n)); end
    if (n == 0) samples = []; bestnorm = 0; samples_sequence = [];
    elseif (n < m) samples = 1:n; bestnorm = mcs(b(:,cols)); samples_sequence = ones(n,1);
    else
        samples = zeros(m, 1);
        flag = zeros(n, 1);
        samples_sequence = zeros(n, m);
        bestnorm = zeros(m, 1);
        for sam = 1 : m
            if (sam == 1)
                if (seed == 0) seed = randi(n); end
                samples(sam) = seed;
                flag(seed) = 1;
                samples_sequence(:,sam) = flag;
                bestnorm(sam) = 1;
                if (debug) fprintf('(%d,B,%.1f):seed %d, norm=1\n', m, rk, seed); end
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
%             snorms = sort(unique(norms));
%             upper = snorms(min(length(snorms), max(1,floor(rk*n))));
%             pos = find(norms==upper);
%             if (random ~= 0) choose = randi(min(length(pos), floor(rk*n)));end
%             if (length(pos) < choose) choose = length(pos); end
%             choose = randi(length(pos));
            
            [~,pos] = sort(norms);
            choose = randi([floor(rk*n)-10, floor(rk*n)+10]);
            if (choose<=0) choose = 1;end
            if (choose>length(pos)) choose = length(pos);end
%             choose = randi(min(length(pos), max(1,floor(rk*n))));
            picked = pos(choose);
            bestnorm(sam) = norms(picked);
            samples(sam) = idx(picked);
            flag(idx(picked)) = 1;
            samples_sequence(:, sam) = flag;
            if (debug) fprintf('(%d,B,%.2f):%d-th minial norm=%.2f, %d picked in round%d\n',...
                                m, rk, choose, bestnorm(sam), idx(picked), sam); end
%             if (debug) fprintf('(%d,B,%.1f):%d unique norms, %d-th minial norm=%f, %d picked (%d ties) in round%d\n',...
%                                 m, rk, length(snorms), choose, bestnorm(sam), idx(picked), length(pos), sam); end
        end
    end
end