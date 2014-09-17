% produce *near* sparse vectors where coefficients decay with power law
%

function [v pos] = nsparsev(n, r, p)
    if (nargin < 2) r = 1; end
    if (nargin < 3) p = 1; end
    v = zeros(n, 1);
    pos = randperm(n);
    v(pos(1)) = r;
    for i = 2 : n
        v(pos(i)) = 2*(rand-0.5) * r * i^(-1/p);
    end
end