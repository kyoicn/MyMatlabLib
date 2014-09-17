function s = sparsev(n, k, r)
if (nargin < 3) r = 1; end

if (k >= n) 
    s = r * rand(n, 1);
elseif (k <= 0)
    s = zeors(n, 1);
else
    s = zeros(n, 1);
    randp = randperm(n);
    s(randp(1 : k)) = r * rand(k, 1);
end