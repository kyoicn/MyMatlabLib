%% rotate marix

function re=rotm(s,d)
si=size(s);
n=si(1);
if(d==90)
    for i=1:n
        for j=1:n
            re(i,j)=s(n-j+1,i);
        end
    end
elseif(d==-90)
    for i=1:n
        for j=1:n
            re(i,j)=s(j,n-i+1);
        end
    end
end