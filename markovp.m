%% markov sequence with p

function s = markovp(n,p)
d = p;
s=zeros(1,n);
for i = 1 : n
    if (i == 1)
        s(i) = 3 * normrnd(0, d);
    else
        s(i) = 3 * normrnd(s(i - 1), d);
    end
end
s=s';