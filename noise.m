function r=noise(n,d)
r = zeros(n, 1);
for i=1:1:n
    r(i)=normrnd(0,d);
end