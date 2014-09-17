% Partial Inversion
% G. Chen, A. Divekar, D. Needell, 
% "Guaranteed sparse signal recovery with highly coherent sensing matrices"

function [c, idx] = partinv(y, a, l, iter)
if (nargin < 4) iter = 100; end
[m, n] = size(a);
c = a' * y;
[~, pos] = sort(abs(c), 'descend');
idx = pos(1 : l);
improve = 1;
k = 0;
while (improve && k < iter)
    pidx = idx;
    part = a(:, idx);
    c_idx = inv(part' * part) * part' * y;
    r = y - part * c_idx;
    idxc = setdiff((1:n)', idx);
    part = a(:, idxc);
    c_idxc = part' * r;
    [~, pos] = sort(abs(c_idxc), 'descend'); %(?)
    idx = pos(1 : l);
    if (length(intersect(idx, pidx)) == length(union(idx, pidx))) improve = 0; end
    k = k + 1;
end
c = zeros(n, 1);
c(idx) = c_idx;
end