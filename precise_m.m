function [recall precision] = precise_m(s1, s2, thr)
[d1 d2] = size(s1);
s1c = abovethr(s1, thr);
s2c = abovethr(s2, thr);
c = 0;
for i = 1:d1
    for j = 1:d2
        if (s2(i, j) > thr && s1(i, j) > thr)
            c = c + 1;
        end
    end
end
recall = c / s1c;
precision = c / s2c;