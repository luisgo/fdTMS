function [Aine,bine]=energyconstraints(LL,pererr,ncoils,ub,lb,NA)
stvar=ncoils+NA;
nsamp=round(sqrt(4*pererr));
Aine=zeros([nsamp*ncoils  2*ncoils+NA]);
bine=50*ones([nsamp*ncoils 1]);
for i=1:ncoils
    Xmax=ub(i);
    Xmin=lb(i);
    for j=1:nsamp
        lam=(j-1/2)/nsamp;
        xoi=Xmin*lam+Xmax*(1-lam);
    Aine(j+(i-1)*nsamp,i)=2*LL(i,i)*xoi;
    Aine(j+(i-1)*nsamp,stvar+i)=-1;
    bine(j+(i-1)*nsamp)=LL(i,i)*xoi^2;
    end
end
end