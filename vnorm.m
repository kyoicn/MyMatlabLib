%%vector norm

function r = vnorm(v, l)
% n = length(v);
r = 0;
if (l == 0)
    v(v == 0) = [];
    r = length(v);
elseif (l > 0)
    r = sum(abs(v).^l);
    r = r^(1 / l);
end