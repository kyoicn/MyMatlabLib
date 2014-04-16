function s = polynoise(n,d)
%% polynomial
for i = 1 : n
%     s(i) = i;
%     s(i) = 0.00000007 * (i + 15) * (i - 66) * (i - 170) * (i - 290);
%     s(i) = 0.000002 * (i - 26) * (i - 364) * (i - 92);
    s(i) = 0.0002 * (i - 26) * (i - 45) * (i - 92) + normrnd(0, d);
%     s(i) = 0.002 * (i + 15) * (i - 66);
%     s(i) = (i/n)^35-(i/n)^3+5;
end
s=s';