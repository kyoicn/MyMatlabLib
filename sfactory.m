%% signal factory

function s = sfactory(n)
%% polynomial
for i = 1 : n
%     s(i) = i;
%     s(i) = 0.00000007 * (i + 15) * (i - 66) * (i - 170) * (i - 290);
%     s(i) = 0.000002 * (i - 26) * (i - 364) * (i - 92);
%     s(i) = 0.000002 * (i - 26) * (i - 364) * (i - 92) + normrnd(0, 0.1);
%     s(i) = 0.002 * (i + 15) * (i - 66);
%     s(i) = (i/n)^35-(i/n)^3+5;
end

%% markov sequence
d = 1;
for i = 1 : n
    if (i == 1)
        s(i) = 1 * normrnd(0, d);
    else
        s(i) = 1 * normrnd(s(i - 1), d);
    end
end

%% vectorized MRF
% s = vectorize(mrf(n, n, 500000), n, n);