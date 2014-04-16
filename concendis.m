function c=concendis(s)
n = length(s);
p = vnorm(s,2)^2;
c = zeros(1, length(s));
acc = 0;

[t,b] = dct4(s);
t_s = sort(abs(t), 1, 'descend');

for i = 1 : length(s)
    acc = acc + t_s(i)^2;
    c(i) = acc/p;
end