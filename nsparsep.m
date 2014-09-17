% check power law parameter of a *near* sparse vector

function [r, p] = nsparsep(v)
n = length(v);
sv = sort(abs(v),'descend');
r = sv(1);
p = 0;
for i = 2 : n
    pp = log(i)/(log(r)-log(sv(i)));
    if (pp > p) p = pp; end
end
end