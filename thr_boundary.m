%% Boundary detection in CSF_thr
% b: a list of indices that form boundary

function [c, b] = thr_boundary(s, thr)
c = 0;
n = length(s);
flag = 0;
for i = 1 : n
    if (s(i) > thr && flag == 0)
        flag = 1;
        c = c + 1;
        b(c) = i;
    elseif (s(i) <= thr && flag == 1)
        flag = 0;
        c = c + 1;
        b(c) = i;
    end
end
if (flag == 1)
    c = c + 1;
    b(c) = n;
end