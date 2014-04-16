function r = abovethr(s, t)
[d1 d2] = size(s);
r = 0;
for i = 1:d1
    for j = 1:d2
        if (s(i, j) > t)
            r = r + 1;
        end
    end
end