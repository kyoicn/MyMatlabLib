%% Determine the KEY POINTS
% s: source signal
% m: number of key points
% method: recovery algorithm

function k = keypoints(s, m, method)
n = length(s);
k = zeros(1, m);
if (strcmp(method, 'GENERAL'))
    [nb, b] = maxminpoints(s);
    if (nb == m)
        k = b;
    elseif(nb > m) % randomly pick up m out of nb
        tempp = randperm(nb);
        k = b(tempp(1:m));
    elseif (nb < m) % we need another (m - nb) points
        add = zeros(1, m - nb);
        if (nb - 1 >= m - nb)
            gap = b(2 : nb) - b(1 : nb - 1);
            dim = size(gap);
            if (dim(1) == 1) sdim = 2;
            elseif (dim(2) == 1) sdim = 1;
            end
            [sortg, sorti] = sort(gap, sdim, 'descend');
            add = floor(0.5 * (b(sorti(1:m-nb)+1) + b(sorti(1:m-nb))));
        else
            randp = randperm(n);
            p = 0;
            for i = 1 : n
                flag = 0;
                for j = 1 : nb
                    if (randp(i) == b(j))
                        flag = 1;
                        break;
                    end
                end
                if (flag == 0)
                    p = p + 1;
                    add(p) = randp(i);
                end
                if (p == m - nb)
                    break;
                end
            end
        end
        k(1 : nb) = b;
        k(nb + 1 : m) = add;
    end
end