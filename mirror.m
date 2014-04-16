%% mirror a matrix

function re=mirror(s)
si=size(s);
n=si(1);
for i=1:n
    for j=1:n
        re(i,j)=s(n-i+1,j);
    end
end