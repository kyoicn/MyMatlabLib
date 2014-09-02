%% mse
function r=my_mse(s1,s2)
n=length(s1);
r=norm(s1-s2,2)^2/n;