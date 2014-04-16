%% MSESR
function r=msesr(s1,s2)
r=mse(s1,s2)/vnorm(s2,2)^2;