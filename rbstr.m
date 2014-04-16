%% produce a random binary string according to a probability p

function s =rbstr(n, p)
s = zeros(n, 1);
m = floor(n * p);
rp = randperm(n);
s(rp(1:m)) = 1;
end