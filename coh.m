%% matrix coherence

function r=coh(m)
cc=m'*m;
m=min(cc);
mm=min(m);
i=0;j=0;
while(i==j)
    ma=max(cc);
    mma=max(ma);
    if(mma==mm) break; end
    [i,j]=find(cc==mma);
    if(i==j) cc(i,j)=mm; end
end
r = mma;