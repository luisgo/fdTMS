function [Anew,jacnew,slacknew,reflnew,tindnew]...
    =generatenewmesh(Aloo,jacloo,slackloo,refl,tind,A,x,ncoils,range,range2)
%%
close all
for i=1:4
load(strcat('barisamp',num2str(i),'.mat'));
ploo{i}=p;
t2ploo{i}=t2p;
end
npts=numel(A(1,1,:));
incval=zeros([npts 1]);
nptnew=0;
for i=1:npts
    incval(i)=(((A(1,:,i)*x(1:ncoils)).^2+(A(2,:,i)*x(1:ncoils)).^2)>=range(1)^2).*...
(((A(1,:,i)*x(1:ncoils)).^2+(A(2,:,i)*x(1:ncoils)).^2)<=range(2)^2);
if incval(i)==1
nptnew=nptnew+4;
else
nptnew=nptnew+1;
end
end
Anew=zeros([2 ncoils nptnew]);
reflnew=zeros([nptnew 1]);
jacnew=zeros([nptnew 1]);
slacknew=zeros([nptnew 1]);
tindnew=zeros([nptnew 1]);
t2p=zeros([nptnew 3]);
etot=zeros([nptnew 1]);
ct=1;
for i=1:npts
if incval(i)==0
tindnew(ct)=tind(i);
reflnew(ct)=refl(i);
Anew(:,:,ct)=Aloo{reflnew(ct)}(1:2,:,tindnew(ct));
jacnew(ct)=jacloo{reflnew(ct)}(tindnew(ct));
t2p(ct,:)=t2ploo{reflnew(ct)}(tindnew(ct),:);
ct=ct+1;
else
for j=1:4
tindnew(ct)=4*(tind(i)-1)+j;
reflnew(ct)=refl(i)+1;    
Anew(:,:,ct)=Aloo{reflnew(ct)}(1:2,:,tindnew(ct));
jacnew(ct)=jacloo{reflnew(ct)}(tindnew(ct));
t2p(ct,:)=t2ploo{reflnew(ct)}(tindnew(ct),:);
ct=ct+1;
end

end
end
for i=1:nptnew
    etot(i)=(Anew(1,:,i)*x(1:ncoils)).^2+(Anew(2,:,i)*x(1:ncoils)).^2;
    slacknew(i)=slackloo{reflnew(i)}(tindnew(i))...
        .*(etot(i)>=(range2(1))^2)+...
        (etot(i)>=(range2(2))^2);
end

save(strcat('res',num2str(max(reflnew(:))),'.mat'));
clf
for i=1:max(reflnew(:))
subplot(1,3,1),
    hold on
trisurf(t2p(reflnew==i,:),ploo{i}(:,1),ploo{i}(:,2),ploo{i}(:,3),slacknew(reflnew==i));
axis equal
subplot(1,3,2),
    hold on
trisurf(t2p(reflnew==i,:),ploo{i}(:,1),ploo{i}(:,2),ploo{i}(:,3),sqrt(etot(reflnew==i)),'edgealpha',0);
colorbar
axis equal
subplot(1,3,3),
    hold on
trisurf(t2p(reflnew==i,:),ploo{i}(:,1),ploo{i}(:,2),ploo{i}(:,3),double(sqrt(etot(reflnew==i))>=50*10^-4));
axis equal
end
pause(5)