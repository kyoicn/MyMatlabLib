%% ANTI HILBERT CURVE
% 1D --> 2D

function r=antihc(s)
n=length(s);
l=round(sqrt(n));
m=round(l/2);
r=zeros(l,l);
if (n==4)
    r(1,1)=s(2);
    r(1,2)=s(3);
    r(2,1)=s(1);
    r(2,2)=s(4);
else
    r(1:m,1:m)=antihc(s(m^2+1:2*m^2));
    r(1:m,m+1:l)=antihc(s(2*m^2+1:3*m^2));
    r(m+1:l,1:m)=mirror(rotm(antihc(s(1:m^2)),90));
    r(m+1:l,m+1:l)=mirror(rotm(antihc(s(3*m^2+1:n)),-90));
end