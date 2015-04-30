% Cost-aware sampling
% sample_coor = ca_sample(m, cost, alp = 1)

function samples = ca_sample2(m, cost, alp)
    if (nargin < 3) alp = 1; end
    dim = size(cost);
    n = prod(dim);
    if (n < m) samples = ones(dim);
    elseif (n == 0) samples = [];
    else
        prob = cost(:).^(-alp)/sum(cost(:).^(-alp)) * m;
        prob(prob > 1) = 1;
        prob(prob < 0) = 0;
        samples = zeros(dim);
        samples(binornd(1, prob, dim) == 1) = 1;
    end
end