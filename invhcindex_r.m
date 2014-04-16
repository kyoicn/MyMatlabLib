%% convert hcindex into (x,y) in d*d matrix

function [x,y] = invhcindex(i, d)
n = d * d;
if (d == 1)
    x = 1;
    y = 1;
else
    if (i <= n / 4)
        [tx, ty] = invhcindex(i, d/2);
        % u/d mirror, rot l
        x = d + 1 - ty;
        y = d/2 + 1 - tx;
    elseif (i > n / 4 && i <= n / 2)
        [tx, ty] = invhcindex(i - n/4, d/2);
        % same
        x = tx;
        y = ty;
    elseif (i > n / 2 && i <= 3/4 * n)
        [tx, ty] = invhcindex(i - n/2, d/2);
        % same
        x = tx;
        y = d/2 + ty;
    elseif (i > 3/4 * n && i <= n)
        [tx, ty] = invhcindex(i - 3/4 * n, d/2);
        % u/d mirror, rot r
        x = ty + d/2;
        y = tx + d/2;
    end
end