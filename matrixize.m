%% un-vectorize

function c = matrixize(s, n1, n2)
c = zeros(n1, n2);
n = length(s);
if (n == n1 * n2)
    for i = 1 : n1
        if (mod(i, 2) == 0)
            c(i, 1 : n2) = s((i * n2) : -1 : ((i - 1) * n2 + 1));
        else
            c(i, 1 : n2) = s(((i - 1) * n2 + 1) : (i * n2));
        end
    end
end