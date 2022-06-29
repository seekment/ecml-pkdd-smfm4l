function [p,Wt] = getOptimalWhWeights(Ks,beta)
[r,c]=size(Ks);
Kss = Ks(:);
[sorted, index] = sort(Kss, 'ascend');
index_zero = find(Kss==0);
index2 = setdiff(index, index_zero,'stable');
newKss = Kss(index2);

K=length(newKss);
p=K;
bfind=1;
miu=0;
while p>0 && bfind
    miu=(sum(newKss(1:p))+2*beta)/p;
    if (miu-newKss(p))>0
        bfind=0;
    else
        p=p-1;
    end
end
% fprintf('miu=%d\n',miu);
newWt=zeros(K,1);
for ii=1:p
    newWt(ii)=(miu-newKss(ii))/(2*beta);
end

Wt=zeros(r,c);
Wt(index2)=newWt;
end
