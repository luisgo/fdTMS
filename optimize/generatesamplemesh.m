function [reflnew,tindnew,pcnew,t2p,ploo]=generatesamplemesh(refl,tind,Emag,range)
%%
close all
for i=1:4
load(strcat('barisamp',num2str(i),'.mat'));
ploo{i}=p;
t2ploo{i}=t2p;
pcen{i}=(p(t2p(:,1),:)+p(t2p(:,2),:)+p(t2p(:,3),:))/3;
end
npts=numel(Emag);
incval=zeros([npts 1]);
nptnew=0;
for i=1:npts
    incval(i)=(Emag(i)<=range(2)).*(Emag(i)>=range(1));
if incval(i)==1
nptnew=nptnew+4;
else
nptnew=nptnew+1;
end
end
reflnew=zeros([nptnew 1]);
tindnew=zeros([nptnew 1]);
t2p=zeros([nptnew 3]);
pcnew=zeros([nptnew 3]);
ct=1;
for i=1:npts
if incval(i)==0
tindnew(ct)=tind(i);
reflnew(ct)=refl(i);
t2p(ct,:)=t2ploo{reflnew(ct)}(tindnew(ct),:);
pcnew(ct,:)=pcen{reflnew(ct)}(tindnew(ct),:);
ct=ct+1;
else
for j=1:4
tindnew(ct)=4*(tind(i)-1)+j;
reflnew(ct)=refl(i)+1;    
t2p(ct,:)=t2ploo{reflnew(ct)}(tindnew(ct),:);
pcnew(ct,:)=pcen{reflnew(ct)}(tindnew(ct),:);
ct=ct+1;
end

end
end
end
