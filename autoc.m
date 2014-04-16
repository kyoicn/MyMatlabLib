%% auto-corr

function r=autoc(s,a)
n=length(s);
s1=zeros(1,n);
for i=1:n
    j=i+a;
    if(j>n)
        j=j-n;
    end
    s1(i)=s(j);
end
re=corrcoef(s,s1);
r=re(1,2);