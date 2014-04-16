%% uniformity of cost distribution

function u = uniformity_cut(cost, line)
    si = size(cost);
    sum1 = 0;
    sum2 = 0;
    count1 = 0;
    
    v = 0;
    h = 0;
    slope = 0;
    
    if (line(1) == line(3))
        % vertical
        v = 1;
    elseif (line(2) == line(4))
        % horizontal
        h = 1;
    else
        slope = (line(4) - line(2)) / (line(3) - line(1));
    end
    for i = 1 : si(1)
        for j = 1 : si(2)
            if (v == 1 && line(4) >= line(2) || h == 1 && line(3) >= line(1) || j >= line(2) + (i - line(1)) * slope)
                sum1 = sum1 + cost(i, j);
                count1 = count1 + 1;
            else
                sum2 = sum2 + cost(i, j);
            end
        end
    end
    if (sum1 * sum2 == 0)
        u = 1;
    else
        u = sum1 / sum2 * (si(1) * si(2) - count1) / count1;
    end
    if (u < 1)
        u = 1 / u;
    end
end