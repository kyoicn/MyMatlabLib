%% convert (x,y) into index in hc d*d

function re=hcindex(d, x, y)
if (d == 1)
    re = 1;
else
    if (x > d/2 && y <= d/2) %I
        x = x - d/2;        
        re = hcindex(d/2, d/2 - y + 1, d/2 - x + 1);
    elseif (x <= d/2 && y <= d/2) %II
        re = d^2 / 4 + hcindex(d/2, x, y);
    elseif (x <= d/2 && y > d/2) %III
        y = y - d/2;
        re = d^2 / 2 + hcindex(d/2, x, y);
    else %IV
        x = x - d/2;
        y = y - d/2;
        re = 3 * d^2 / 4 + hcindex(d/2, y, x);
    end
end