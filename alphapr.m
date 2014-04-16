%% alpha-aware probability generator

function pr = alphapr(cost, alp)
    pr = cost.^(-alp) / sum(sum(cost.^(-alp)));
end