%% Matrix completion

function x = matrixcompletion(measurement, indicator, rank, lamda, t)
    si = size(measurement);
    m = si(1);
    n = si(2);
    L = rand(m, rank);
    sI = sqrt(lamda) * eye(rank);
    
    bestV = 1e99;
    bestL = L;
    bestR = rand(rank, n);
    for i = 1 : t
        tP = [L; sI];
        tQ = [measurement; zeros(rank, n)];
        R = (tP' * tP) \ tP' * tQ;
        tP = [R'; sI];
        tQ = [measurement'; zeros(rank, m)];
        L = (tP' * tP) \ tP' * tQ;
        L = L';
        v = norm(indicator .* (L * R) - measurement, 'fro')^2;
        v = v + lamda * (norm(L, 'fro')^2 + norm(R', 'fro')^2);
        if (v < bestV)
            bestV = v;
            bestL = L;
            bestR = R;
        end
    end
    x = bestL * bestR;
end