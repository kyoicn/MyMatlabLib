function r = errpluscost(err, cost)
    l = 0.5;
    r = l * err + (1 - l) * cost;
end