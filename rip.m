function [re, deltas] = rip(a, k)
re = 0;
si = size(a);
if (si(1) >= si(2))
    n = si(1);
    a = a';
else
    n = si(2);
end

iter = 100;
deltas = zeros(1, iter);
for i = 1 : iter
    sv = zeros(1, n);
    randp = randperm(n);
    for j = 1 : k
        sv(randp(j)) = rand();
    end
    delta = abs(vnorm(a * sv', 2)^2 / vnorm(sv, 2)^2 - 1);
    if re < delta
        re = delta;
    end
    deltas(i) = delta;
end