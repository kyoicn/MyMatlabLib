%% uniformity of cost distribution

function [u_avg start] = uniformity_win(cost, win)
    si = size(cost);
    n = si(1) * si(2);
    rest = n - win^2;
    costsum = sum(sum(cost));
    start = [1 1];
    u_avg = 1;
    for i = 1 : si(1) - win + 1
        for j = 1 : si(2) - win + 1
            block = cost(i : i + win - 1, j : j + win - 1);
            blocksum = sum(sum(block));
            restavg = (costsum - blocksum) / rest;
            ratio_avg = blocksum / win^2 / restavg;
            if (ratio_avg < 1)
                ratio_avg = 1 / ratio_avg;
            end
            if (ratio_avg > u_avg)
                u_avg = ratio_avg;
                start = [i j];
            end
        end
    end
end