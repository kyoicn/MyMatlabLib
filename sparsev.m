function s = sparsev(n, k)
if (k >= n) 
    s = 100 * rand(n, 1);
elseif (k <= 0)
    s = zeors(n, 1);
else
    s = zeros(n, 1);
    randp = randperm(n);
    s(randp(1 : k)) = 100 * rand(k, 1);
end