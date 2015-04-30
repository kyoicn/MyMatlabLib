% calculate (absolute) column sum

function [max_cs origin_cs cs_array] = mcs(b)
    s = sum(b, 1);
    cs_array = abs(s);
    [max_cs, pos] = max(cs_array);
    origin_cs = s(pos);
end