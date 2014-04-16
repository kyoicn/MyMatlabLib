%% sparsity level of a vector

function [sc, in] = sparsec(v, alp)
    if nargin < 2
        alp = 0.95;
    end
    n = length(v);
    n2 = norm(v)^2;
    [sv, idx] = sort(abs(v), 'descend');
    ac = 0;
    for i = 1 : n
        ac = ac + sv(i)^2;
        if ac / n2 >= alp
            break;
        end
    end
    sc = i;
    in = idx(1 : sc);
end