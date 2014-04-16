%% hamming weight

function hw = hweight(s);
si = size(s);
di = 0;
if (si(2) == 1)
    s = s';
    di = si(1);
else
    di = si(2);
end
hw = di * pdist([s; zeros(1, di)], 'hamming');