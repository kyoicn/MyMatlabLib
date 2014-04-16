%GMP
function x_hat = gmp(y, a, lower, upper)
Size = size(a);
m = Size(1);
n = Size(2);

count = n;
x_hat = zeros(n, 1);
N = 1 : n;
r_y = y;
z = zeros(n, 1);
while (1)
    norms = [];
    pos = 0;
    for i = 1 : n
        if (N(i) ~= -1)
            pos = pos + 1;
            tmp = (a(:, i)' * r_y) / (a(:, i)' * a(:, i));
            if (tmp <= upper & tmp >= lower)
                v = zeros(n, 1);
                v(i) = tmp;
                norms(1:3, pos) = [norm(r_y - a * v, 2); i; tmp];
            else
                lv = zeros(n, 1);
                uv = zeros(n, 1);
                lv(i) = lower;
                uv(i) = upper;
                if (norm(r_y - a * lv, 2) <= norm(r_y - a * uv, 2))
                    norms(1:3, pos) = [norm(r_y - a * lv, 2); i; lower];
                else
                    norms(1:3, pos) = [norm(r_y - a * uv, 2); i; upper];
                end
            end
        end
    end
    if (count == 0)
        break;
    end
    [min_norm, index] = min(norms(1, :));
    i = norms(2, index);
    z(i) = norms(3, index);
    if (z(i) == 0)
        break;
    end
    N(i) = -1;
    count = count - 1;
    x_hat(i) = z(i);
    r_y = r_y - a * z;
    z = zeros(n, 1);
end