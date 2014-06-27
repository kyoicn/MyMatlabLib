% Sparsity level of a vector,
% or energy concentrated on the major coefficients

function [sc, in] = sparsec(v, alp)
    if nargin < 2
        alp = 0.95;
    end
    n = length(v);
    [sv, idx] = sort(abs(v), 'descend');
    n2 = norm(v)^2;
    
    if (alp <= 1)
        ac = 0;
        for i = 1 : n
            ac = ac + sv(i)^2;
            if ac / n2 >= alp
                break;
            end
        end
        sc = i;
        in = idx(1 : sc);
        
    elseif (alp > 1)
        alp = floor(alp);
        partial = v(idx(1:alp));
        sc = norm(partial)^2 / n2;
        in = idx(1 : alp);
    end
        
end