function [Aine,bine,NA,slackid]=peakEfieldconstraintlinear(ncoils,Nang,slack,A,condi,stimth,maxE)
nslack=numel(slack);
NA=nnz(slack==1);
thi=0:2*pi/Nang:2*pi-2*pi/Nang;%discretize angles in the circle
Aine=zeros([Nang*nslack  ncoils+NA]);
bine=stimth*ones([Nang*nslack 1]);
%write field constraints
ct=1;
slackid=zeros([Nang*NA 1]);
for i=1:nslack
    th=thi+rand(1)*2*pi/Nang;%shift it so points dont line up
for j=1:Nang
    Aine(j+(i-1)*Nang,1:ncoils)=(A(1,:,i)*sin(th(j))+A(2,:,i)*cos(th(j)))*condi;
        if slack(i)==1
    Aine(j+(i-1)*Nang,ncoils+ct)=-(maxE-stimth);
    slackid(j+(ct-1)*Nang)=j+(i-1)*Nang;
        end
if slack(i)==2
bine(j+(i-1)*Nang)=maxE;
end
end
if slack(i)==1
ct=ct+1;
end
end
Aine=sparse(Aine);
end