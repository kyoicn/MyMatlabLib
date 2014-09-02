% Minimal norm sampling: with cost and reference column sum bounds
% samples = minnorm_sample(m, base)
%
% m: number of samples
% b: base
% id: identity column index
% rk:

function [samples, bestnorm, samples_sequence] = minnorm_sample_cost_ref(m, b, id, cost, ref, seed, random, debug)
    % so far only deal with 2D
    if (nargin < 6) seed = 0; random = 0; debug = 0;
    elseif (nargin < 7) random = 0; debug = 0;
    elseif (nargin < 8) debug = 0;
    end
    dim = size(b);
    n = dim(2);
    sumcost = sum(cost);
    acost = 0;
    cols = 1:n;
    cols(id) = [];
    if (debug) normf = figure; end
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
                if (debug) fprintf('(%d,B=%.2f):seed %d, norm=%.2f, cost=%.1f%%(%.1f%%)\n',...
                        m, ref(sam), seed, bestnorm(sam), acost/sumcost*100, (acost/sumcost-sam/n)*100); end
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
            % col sum map
            allnorm = nan(1, n);
            allnorm(idx) = norms;
            if (debug && sam < n)
                figure(normf);
                hold off;
                plot(allnorm);
                hold on;
                plot(samplex, min(norms)*ones(1, length(samplex)), 'o', 'MarkerSize', 5, 'MarkerFaceColor', 'k');
                plot(1:n, ref(sam)*ones(1, n), 'r', 'LineWidth', 2);
                xlim([0 n+1]);
                ylim([min(norms)-1e-10 max(norms)+1e-10]);
            end
            %normmap = antihc(allnorm);
            
%             [~, pos] = sort(norms);
            pos = find(norms<=ref(sam));
            wflag = 0;
            if (isempty(pos))
                wflag = 1;
                [~, choose] = min(norms);
                picked = choose;
                c = cost(idx(choose));
            else
                [c, choose] = min(cost(idx(pos)));
                picked = pos(choose);
            end
            bestnorm(sam) = norms(picked);
            samples(sam) = idx(picked);
            acost = acost + c;
            flag(idx(picked)) = 1;
            samples_sequence(:, sam) = flag;
            if (debug) 
                fprintf('(%d,B=%.2f): norm=%f, cost=%.1f%%(%.1f%%), %d picked in round%d',...
                        m, ref(sam), bestnorm(sam), acost/sumcost*100, (acost/sumcost-sam/n)*100, idx(picked), sam);
                if (wflag) fprintf('(0)'); end
                fprintf('\n');
            end
                            
        end
    end
end