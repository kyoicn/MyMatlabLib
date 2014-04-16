%% markov sequence with d

function s = markovd(n,d)
s=zeros(1,n);
for i = 1 : n
    if (i == 1)
        s(i) = 2 * normrnd(0, d);
    else
        s(i) = 1 * normrnd(s(i - 1), d);
    end
end
s=s';