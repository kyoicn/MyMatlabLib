%% accuracy score

function acc = acc100(s, r)
    error = norm(s - r, 2);
    if (error > norm(s, 2)) 
        acc = 0;
    else
        if (norm(s, 2) == 0)
            acc = 1;
        else
            acc = 1 - (error / norm(s, 2))^2;
        end
    end
end