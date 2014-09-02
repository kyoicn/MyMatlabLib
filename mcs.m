% calculate (absolute) column sum

function [max_cs cs] = mcs(b)
    cs = abs(sum(b,1));
    max_cs = max(cs);
end