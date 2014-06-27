% Cost-aware sampling
% sample_coor = ca_sample(m, cost, alp = 1)

function samples = ca_sample(m, cost, alp)
    if (nargin < 3) alp = 1; end
    dim = size(cost);
    n = prod(dim);
    if (n < m) samples = ones(dim);
    elseif (n == 0) samples = [];
    else
        samples = zeros(dim);
        samples(randsamplewtr(n, m, cost(:).^(-alp)/sum(cost(:).^(-alp)))) = 1;
    end
end