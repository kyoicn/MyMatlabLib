%% convert hcindex into (x,y) in d*d matrix

function [x,y] = invhcindex(i, d)
    xstack = zeros(1, d);
    ystack = zeros(1, d);
    tstack = zeros(1, d);
    pt = 0;

    while (d > 0)
        n = d * d;
        pt = pt + 1;
        if (d == 1)
            xstack(pt) = 1;
            ystack(pt) = 1;
        else
            if (i <= n / 4)
                xstack(pt) = d + 1;
                ystack(pt) = d / 2 + 1;
                tstack(pt) = 1;
            elseif (i > n / 4 && i <= n / 2)
                xstack(pt) = 0;
                ystack(pt) = 0;
                tstack(pt) = 2;
                i = i - n / 4;
            elseif (i > n / 2 && i <= 3/4 * n)
                xstack(pt) = 0;
                ystack(pt) = d / 2;
                tstack(pt) = 3;
                i = i - n / 2;
            elseif (i > 3/4 * n && i <= n)
                xstack(pt) = d / 2;
                ystack(pt) = d / 2;
                tstack(pt) = 4;
                i = i - 3 / 4 * n;
            end
        end
        d = floor(d / 2);
    end

    pt = pt - 1;
    while (pt > 0)
        if (tstack(pt) == 1)
            xstack(pt) = xstack(pt) - ystack(pt + 1);
            ystack(pt) = ystack(pt) - xstack(pt + 1);
        elseif (tstack(pt) == 2)
            xstack(pt) = xstack(pt + 1);
            ystack(pt) = ystack(pt + 1);
        elseif (tstack(pt) == 3)
            xstack(pt) = xstack(pt + 1);
            ystack(pt) = ystack(pt) + ystack(pt + 1);
        elseif (tstack(pt) == 4)
            xstack(pt) = xstack(pt) + ystack(pt + 1);
            ystack(pt) = ystack(pt) + xstack(pt + 1);
        end
        pt = pt - 1;
    end

    x = xstack(1);
    y = ystack(1);
end