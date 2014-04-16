%% mse
function r=mse(s1,s2)
n=length(s1);
r=vnorm(s1-s2,2)^2/n;