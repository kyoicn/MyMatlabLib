function r = l2normratio(s, r)
    r = norm(s, 2) / norm(s - r, 2);
end