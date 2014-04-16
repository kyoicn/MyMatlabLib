function re = err01(s, r)
    if (snr(s, r) >= 35)
        re = 0;
    else
        re = 1;
    end
end