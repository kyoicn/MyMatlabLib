%% Hilber curve

function re=hilbertcurve(s);
si=size(s);
n=si(1);
if(n==2)
    re=[s(2,1),s(1,1),s(1,2),s(2,2)];
else
    re=[hilbertcurve(rotm(mirror(s((n/2+1):n,1:(n/2))),-90)),hilbertcurve(s(1:(n/2),1:(n/2))),hilbertcurve(s(1:(n/2),(n/2+1):n)),hilbertcurve(rotm(mirror(s((n/2+1):n,(n/2+1):n)),90))];
end
