function [adj count] = thr_adj(s2,thr)
[d1, d2] = size(s2);
adj = zeros(d1, d2);
count = 0;
for i = 1:d1
    for j = 1:d2
        if (s2(i, j) > thr)
            count = count + 1;
            adj(i, j) = s2(i, j);
        else
            adj(i, j) = thr;
        end
    end
end