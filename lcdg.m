%cdg luo

function [re, mse] = lcdg(s, n, m, k)
% n = 100;
% m = 50;
% k = 10;
% 
% originData = sfactory(n);

for i = 1 : m
    for j = 1 : n
        phi(i, j) = rand();
    end
end

[t, b] = dct4(s);

a = phi * b^(-1);
% r = omp(phi * s', a, n, k);
r=bp(phi*s',a,n);
re = b^(-1) * r';

% figure;
% plot(1 : n, s, 'k-', 'LineWidth', 2);
% hold on;
% plot(1 : n, re, 'r-.', 'LineWidth', 3);
% 
% figure;
% bar(1 : n, t);
% hold on;
% plot(1 : n, r, 'r*');

mse = 0;
for i = 1 : n
    mse = mse + (re(i) - s(i))^2;
end
mse = mse / n;