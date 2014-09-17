% calculate (absolute) column sum

function [max_cs cs] = amcs(b)
    cs = abs(sum(b,1))/size(b,1);
    max_cs = max(cs);
end