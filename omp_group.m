%OMP Recovery
%
% stop: stop condition

function co = omp_group(y, A, n, k, stop)

si = size(A);
m = si(1);
si = size(y);
count = si(2);
co = zeros(n, count);
residual = y;
pos_array = zeros(k, count);
AProduct = A' * A;

diagy = zeros(m * count, count);
for i = 1 : count
    diagy((i - 1) * m + 1 : i * m, i) = y(:, i);
end

% Aug_t = zeros(m, k);
% mask = zeros(n, count);
% oneton = 1 : count;

for times = 1 : k
    product = abs(A' * residual);
    [~, pos] = max(product);
    pos_array(times, :) = pos;
%     mask([pos' oneton']) = ones(1, count);
    toinv = zeros(times * count);
    diagpart = zeros(times * count, m * count);
    concatpart = zeros(m, times * count);
    for i = 1 : count
        toinv((i - 1) * times + 1: i * times, (i - 1) * times + 1 : i * times) = AProduct(pos_array(1 : times, i), pos_array(1 : times, i));
        diagpart((i - 1) * times + 1 : i * times, (i - 1) * m + 1 : i * m) = A(:, pos_array(1 : times, i))';
        concatpart(:, (i - 1) * times + 1 : i * times) = A(:, pos_array(1 : times, i));
    end
    diagparty = toinv \ diagpart * diagy;
    residual = y - concatpart * diagparty;
    
    if (stop < 0)
        continue;
    elseif (all(norm(residual) <= stop))
        break;
    end
end

if (times < k)
    % remove zeros
    pos_array = pos_array(1 : times);
end
    
for i = 1 : count
    co(pos_array(:, i), i) = diagparty((i - 1) * times + 1 : i * times, i);
end
