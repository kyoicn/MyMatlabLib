%% matrix coherence

function r=coh(m)
dim = size(m);
cc=abs(m'*m);
diagelements = diag(cc);
p1 = repmat(sqrt(diagelements), 1, dim(2));
c = cc.*(ones(dim(2))-diag(ones(1,dim(2))))./p1./p1';
r = max(c(:));
% r = max(abs(m(:)));
% m=min(cc);
% mm=min(m);
% i=0;j=0;
% while(i==j)
%     ma=max(cc);
%     mma=max(ma);
%     if(mma==mm) break; end
%     [i,j]=find(cc==mma);
%     if(i==j) cc(i,j)=mm; end
% end
% r = mma;