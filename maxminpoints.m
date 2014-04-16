function [nb, b] = maxminpoints(s)
n = length(s);
nb = 0;
flag = 0;
for i = 1:n
    if (i == 1)
        nb = nb + 1;
        b(nb) = i;
        if (n >= 2)
            if (s(2) >= s(1)) flag = 0;
            else flag = 1;
            end
        end
    elseif (i == n)
        nb = nb + 1;
        b(nb) = i;
    else
        if (s(i + 1) > s(i) && flag == 1)
            flag = 0;
            nb = nb + 1;
            b(nb) = i;
        elseif (s(i + 1) < s(i) && flag == 0)
            flag = 1;
            nb = nb + 1;
            b(nb) = i;
        end
    end
end