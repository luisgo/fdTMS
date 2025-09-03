function [t2e,e2t,e2p,bnodes,inodes,Aloop,Arwg,tvol,evol,nhat]=extractmesharrays(t2p,p);
debug=0;
nt=size(t2p);
np=size(p);
edges=[t2p(:,[2,3]);t2p(:,[1,3]);t2p(:,[1,2]);]; %all possible edges with repeticions
node6=[(1:nt(1))';(1:nt(1))'+nt(1);(1:nt(1))'+2*nt(1)];%%indexing for edges
node5=[(1:nt(1)); (1:nt(1)); (1:nt(1)); ]';
edges=sort(edges,2);%allways lead with lowest point
[e2p,~,jx]=unique(edges,'rows','stable'); %find unique edges
t2e(node6)=jx;
t2e=reshape(t2e,[nt(1) 3]);
ne=size(e2p);
e2t=zeros([ne(1) 2]);
for i=1:length(jx) %extracts appropriate triangle
if e2t(jx(i),1)==0
    e2t(jx(i),1)=node5(i);
elseif e2t(jx(i),2)==0
    e2t(jx(i),2)=node5(i);
end
end


volfunc= @(v1,v2) sqrt((v1(:,2).*v2(:,3)-v1(:,3).*v2(:,2)).^2+...
                       (v1(:,3).*v2(:,1)-v1(:,1).*v2(:,3)).^2+...
                       (v1(:,1).*v2(:,2)-v1(:,2).*v2(:,1)).^2)/2;
 tvol=volfunc(p(t2p(:,2),:)-p(t2p(:,1),:),p(t2p(:,3),:)-p(t2p(:,2),:));
 evol=sqrt(sum(p(e2p(:,2),:)-p(e2p(:,1),:).^2,2));
% 
trimesh(t2p,p(:,1),p(:,2),p(:,3),'facealpha',0);
    tcen=zeros([nt(1) 3]);
    nhat=zeros([nt(1) 3]);
for i=1:nt(1)
    nh=cross(p(t2p(i,2),:)-p(t2p(i,1),:),p(t2p(i,3),:)-p(t2p(i,1),:));
    nhat(i,:)=nh/norm(nh);
    hold on
    tcen(i,1)=sum(p(t2p(i,:),1),1)/3;
    tcen(i,2)=sum(p(t2p(i,:),2),1)/3;
    tcen(i,3)=sum(p(t2p(i,:),3),1)/3;
if debug==1
quiver3(tcen(i,1),tcen(i,2),tcen(i,3),nh(1),nh(2),nh(3));
text(tcen(i,1),tcen(i,2),tcen(i,3),num2str(i),'Fontsize',12);
end
end
    ecen=zeros([ne(1) 3]);
for i=1:ne(1)
    ecen(i,1)=sum(p(e2p(i,:),1),1)/2;
    ecen(i,2)=sum(p(e2p(i,:),2),1)/2;
    ecen(i,3)=sum(p(e2p(i,:),3),1)/2;
    if debug==1
    text(ecen(i,1),ecen(i,2),ecen(i,3),num2str(i),'Fontsize',12);
    va=(tcen(e2t(i,1),:)+ecen(i,:))/2;
    text(va(1),va(2),va(3),num2str(1),'Fontsize',12);
    if e2t(i,2)~=0
    va=(tcen(e2t(i,2),:)+ecen(i,:))/2;
    text(va(1),va(2),va(3),num2str(2),'Fontsize',12);
    end
    end
end
npt=size(p);
pcen=zeros([ne(1) 3]);
for i=1:npt(1)
    pcen(i,1)=sum(p(i,1),1);
    pcen(i,2)=sum(p(i,2),1);
    pcen(i,3)=sum(p(i,3),1);
    if debug==1
    text(pcen(i,1),pcen(i,2),pcen(i,3),num2str(i),'Fontsize',24);
    end
end
hold on;
% 
col=ones([2*ne(1) 1]);
row=ones([2*ne(1) 1]);
val=zeros([2*ne(1) 1]);
nnodes = 0;
ientry=0;
np=size(p);
n2p=zeros([np(1) 1]);
nodefound=zeros([ne(1) 2]);
for iloop=1:2
    
for iedget = 1:ne
if (iloop==1 && e2t(iedget,2)==0) || (iloop==2 && e2t(iedget,2)~=0) %iloop one only does boundary points
for inode = 1:2
    
if nodefound(iedget, inode) == 0  % 0 is false
    if inode==1
    na=e2p(iedget,inode);nb=e2p(iedget,2);
    ehat=p(nb,:)-p(na,:);ehat=ehat/norm(ehat);
    elseif inode==2
    na=e2p(iedget,inode);nb=e2p(iedget,1);
    ehat=p(nb,:)-p(na,:);ehat=ehat/norm(ehat);
    end
     nodefound(iedget, inode) = 1;  % 1 is true
	 nnodes=nnodes+1;%its defining a nodal basis
     if iloop==1
     n2p(na)=na;
     end
     t1t2 = e2t(iedget,:); %go from one triangle to the next
     [~,ph]=intersect(t2e(t1t2(1),:),iedget);
     if t1t2(2)==0
         tnext = t1t2(1);
         firstsign=-1;
            dir1=-sign(sum((tcen(t1t2(1),:)-p(t2p(t1t2(1),ph),:))...
                .*cross(nhat(t1t2(1),:),ehat)));
     else
         tnext = t1t2(2);
         firstsign=1;
         
            dir1=sign(sum((tcen(t1t2(1),:)-p(t2p(t1t2(1),ph),:))...
                .*cross(nhat(t1t2(1),:),ehat)));

     end
     enext=iedget;
     ientry =ientry+1; col(ientry)=enext; row(ientry)=na; val(ientry)=dir1*firstsign;  %write matrix entry

loopcomplete = 0;
while (loopcomplete == 0)
    
for i=1:3
if t2e(tnext,i)~=enext && t2p(tnext,i)~=na%anothe edge with na
    enext=t2e(tnext,i);
   if e2p(enext,1)==na;
       edgefoundind=1;
   elseif e2p(enext,2)==na;
      edgefoundind=2;
   end
   break;
end
end 

nodefound(enext,edgefoundind)=1;

 t1pt2p = e2t(enext,:);
if (enext == iedget || t1pt2p(2) == 0)%either loop is done or boundary
    loopcomplete = 1;
end
if enext~=iedget %move to next triangle
    if (t1pt2p(1) == tnext) 
        tnext = t1pt2p(2);
        secondsign =1;
    elseif (t1pt2p(2) == tnext)
        tnext = t1pt2p(1);
         secondsign=-1;
         
    end
     ientry =ientry+1; col(ientry)=enext; row(ientry)=na; val(ientry)=dir1*secondsign;  %write matrix entry 
end
end
end
end
end
end
end
ppp=1:np(1);
bnodes=n2p(n2p~=0);
inodes=ppp(n2p==0);
Aloop=sparse(col,row,val);


col=ones([2*ne(1) 1]);
row=ones([2*ne(1) 1]);
val=zeros([2*ne(1) 1]);
ct=1;
for i=1:ne
[~,ia]=intersect(t2e(e2t(i,1),:),i);
col(ct)=ia+(e2t(i,1)-1)*3;
row(ct)=i;
val(ct)=1;
ct=ct+1;
if e2t(i,2)~=0
[~,ia]=intersect(t2e(e2t(i,2),:),i);
col(ct)=ia+(e2t(i,2)-1)*3;
row(ct)=i;
val(ct)=-1;
ct=ct+1;
end
end
Arwg=sparse(col,row,val);

