% can we trust the sparsity level in the recovery by inputting
% A, and A*x?

function [yes, ratio] = verifyk(y, a, vanishing)
%-----------
% default vanishing = 20
vanishing = abs(vanishing);
%-----------
si = size(a);
m = si(1);
n = si(2);
upperK = floor(m / log2(n));
x = bp(y, a, n);
[~, in] = sort(abs(x), 'descend');
x = x(in);
x2 = zeros(n, 1);
x2(1:n - 1) = x(2:n);
x2(n) = x(1);
ratio = x ./ x2;
actualK = -1;
for i = 1 : n - 1
    if (ratio(i) > vanishing || ratio(i) < -1 * vanishing)
        actualK = i;
        break;
    end
end
if (actualK == -1)
    actualK = n;
end
if (actualK <= upperK)
    yes = 1;
else
    yes = 0;
end