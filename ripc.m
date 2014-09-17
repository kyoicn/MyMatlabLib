function [re, deltas] = ripc(a, k, iter, t_max)
if (nargin < 3) iter = 100; end
if (nargin < 4) t_max = iter; end
iter = round(iter);

si = size(a);
if (si(1) >= si(2))
    n = si(1);
    a = a';
else
    n = si(2);
end

deltas = zeros(1, iter);
%init
count = 1;
sv = sparsev(n, k);
supp = find(sv~=0);
unsupp = find(sv==0);
ndelta = abs(norm(a*sv)^2/norm(sv)^2 - 1);
cdelta = ndelta;
deltas(count) = cdelta;
temp = 0;
nochange = 0;

while (count < iter && nochange < t_max)
    bak_sv = sv;
    bak_supp = supp;
    bak_unsupp = unsupp;
    
    % transit
    ex1 = randi(length(supp));
    ex2 = randi(length(unsupp));
    sv(supp(ex1)) = 0;
    sv(unsupp(ex2)) = rand();
    tmp = supp(ex1);
    supp(ex1) = unsupp(ex2);
    unsupp(ex2) = tmp;
    
    ndelta = abs(norm(a*sv)^2/norm(sv)^2 - 1);
    if (ndelta > cdelta) accept = 1;
    else
        accept = binornd(1, exp((ndelta-cdelta)/temp));
    end
    
    if (accept)
        count = count + 1;
        cdelta = ndelta;
        deltas(count) = cdelta;
        temp = count / t_max;
    else
        sv = bak_sv;
        supp = bak_supp;
        unsupp = bak_unsupp;
        nochange = nochange + 1;
    end
end
re = max(deltas);