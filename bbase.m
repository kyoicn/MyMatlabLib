%% bbase

function b = bbase(n)
b = zeros(n,n);
for i = 1 : n
    flag = 1;
    for j = 1 : n
        if (flag == 1)
            b(i, j) = 1;
        else
            b(i, j) = 0;
        end
        if (i == 1)
            flag = flag * (-1);
        else
            if (mod(j, i) == 0)
                flag = flag * (-1);
            end
        end
    end
end