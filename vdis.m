function dis = vdis(s, step)
n = length(s);
mi = min(s);
ma = max(s);
dis = zeros(2, floor((ma - mi) / step) + 1);
dis(1, :) = mi : step : ma;
for i = 1 : n
    j = floor((s(i) - mi) / step) + 1;
    dis(2, j) = dis(2, j) + 1;
end
dis(2, :) = dis(2, :) / n;