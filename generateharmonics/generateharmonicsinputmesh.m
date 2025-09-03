function generateharmonicsinputmesh(outputfile,nlayers,meshfile)
debug=2;
load(meshfile)
[pdd,~,t2pp]=unique(t2p(:));
p=p(pdd,:);
t2p=reshape(t2pp,size(t2p));
if debug==2
    clear p t2p;

p(1,:)=[0,0,0];
p(2,:)=[-1,0,0];
p(3,:)=[0,-1,0];
p(4,:)=[-1,-1,0];
p(5,:)=[-1,-2,0];
p(6,:)=[0,-2,0];
t2p=delaunay(p(:,1),p(:,2));

end
if debug>=1;
    disp('if mesh displayed is correct press enter, otherwise debug!')

trisurf(t2p,p(:,1),p(:,2),p(:,3));%plot mesh
pause
close all
end

%% this ensures mesh is in the correct quadrant
%%%processing might need to change for newer meshes
if nnz(p(:,2)>0)>nnz(p(:,2)<0)
p(:,2)=-p(:,2);%force to be in the 4th quadrant
end
if nnz(p(:,1)>0)>nnz(p(:,1)<0)
p(:,1)=-p(:,1);%force to be in the 4th quadrant
end
%make zeros actual zeros
p(abs(p(:,2))<3*10^-5,2)=0;
p(abs(p(:,1))<3*10^-5,1)=0;
%remove unused points
[pind,~,tp]=unique(t2p(:));
p=p(pind,:);
t2p=reshape(tp,size(t2p));
%add a second layer
t2p2=t2p+numel(p(:,3));
if nlayers==2
p2=p;
p2(:,3)=p(:,3)+wireth;
t2p=cat(1,t2p,t2p2);
p=cat(1,p,p2);
end

if debug==1;
    disp('multilayer mesh mesh displayed is correct press enter, otherwise debug!')

trisurf(t2p,p(:,1),p(:,2),p(:,3));%plot mesh
pause
close all
end

%%%%generate that information and nhat information
that1=t2p;
that2=t2p;
nhat=t2p;
for i=1:numel(t2p(:,1))
nh=cross(p(t2p(i,2),:)-p(t2p(i,1),:),p(t2p(i,3),:)-p(t2p(i,1),:));
that1(i,:)=p(t2p(i,2),:)-p(t2p(i,1),:);
that1(i,:)=that1(i,:)/norm(that1(i,:));
that2(i,:)=cross(that1(i,:),nh);
that2(i,:)=that2(i,:)/norm(that2(i,:));
nhat(i,:)=nh;
if nh(3)<0
t2p(i,:)=t2p(i,[2;1;3]);
end
end
if debug==2
nhat
that1
that2
end
%% Generate loop basis sets
%%write symmetry matrix
close all
hold on
np=numel(p(:,1));
p2=zeros([np*4 3]);

nt=numel(t2p(:,1));
t2p2=zeros([nt*4 3]);
s(:,1)=[1;1;-1;-1];
s(:,2)=[1;-1;1;-1];
s(:,3)=[-1;1;-1;1];
col=zeros([4*np 1]);
row=zeros([4*np 1]);
val=zeros([4*np 1]);
ct=1;
for j=1:4
    p2(np*(j-1)+(1:np),1)=s(j,1)*p(:,1);
    p2(np*(j-1)+(1:np),2)=s(j,2)*p(:,2);
    p2(np*(j-1)+(1:np),3)=p(:,3);
    t2p2(nt*(j-1)+(1:nt),:)=t2p+(j-1)*np;
    col(ct:ct+np-1)=np*(j-1)+(1:np);
    row(ct:ct+np-1)=1:np;
    val(ct:ct+np-1)=s(j,3);
    ct=ct+np;
end

%remove redundant points
[psym,ia,ic1] = unique(p2,'rows');
t2psym=ic1(t2p2);
Asym=sparse(col,row,val);

axis equal
Asym=Asym(ia,:);

if debug==2
Asym
trisurf(t2psym,psym(:,1),psym(:,2),psym(:,3));
hold on
for i=1:numel(psym(:,1))
text(psym(i,1),psym(i,2),psym(i,3),num2str(i));
end
pause
end
[t2e,e2t,e2p,bnodes,inodes,Aloop,Arwg,tvol,evol,nhat]=extractmesharrays(t2p,p);
del=2*max(sqrt(2*abs(tvol(:))))
if debug==2
inodes
nhat
Aloop
Arwg
tvol
nhat
end
cctt=1;
keep2=zeros([numel(p(:,1)),1]);
keep2(inodes)=1;
 for i=1:numel(p(:,1))
     if keep2(i)==1
         continue;
     elseif p(i,1)==0 && p(i,2)~=0 && p(i,2)~=min(p(:,2))
     keep2(i)=1;
     end
 cctt=cctt+1;
 end
 

keep=1:numel(p(:,1));
keep=keep(keep2==1);
if debug==2
inodes
nhat
Aloop
Arwg
tvol
nhat
keep
close all
end
if debug>=1
scatter3(p(keep,1),p(keep,2),p(keep,3))
hold on
trisurf(t2p,p(:,1),p(:,2),p(:,3))
cctt=1;
pause
close all
end

%% compute dipole approximations

nt=size(t2p);
% construct from half-RWG to quadrature weight matrices
% nquad=4;
% a=3/5;b=1/5;
% q(:,1)=[a;b;b];
% q(:,2)=[b;a;b];
% q(:,3)=[b;b;a];
% q(:,4)=[1/3;1/3;1/3];
% qwt(1:3)=25/96*2;
% qwt(4)=-9/32*2;
nquad=1;
q(:,1)=[1;1;1]/3;
qwt=1;
nt=nt(1);
Ax=zeros([nquad*nt*4 3*nt]);
Ay=zeros([nquad*nt*4 3*nt]);
Az=zeros([nquad*nt*4 3*nt]);

Apx=zeros([nquad*nt*4 3*nt]);
Apy=zeros([nquad*nt*4 3*nt]);
Apz=zeros([nquad*nt*4 3*nt]);
loc=zeros([nt*nquad*4 3]);
ct=1;
clear s
s(:,1)=[1;1;-1;-1];
s(:,2)=[1;-1;1;-1];
s(:,3)=[-1;-1;1;1];

for ori=1:4
for i=1:nt 
pts=p(t2p(i,:),:);
pts(:,1)=s(ori,1)*pts(:,1);
pts(:,2)=s(ori,2)*pts(:,2);
    for k=1:nquad
  loci=k+(i-1)*nquad+(ori-1)*nquad*nt;
loc(loci,:)=pts(1,:)*q(1,k)+pts(2,:)*q(2,k)+pts(3,:)*q(3,k);
    for j=1:3
Ax(loci,j+(i-1)*3)=s(ori,3)*qwt(k)*(loc(loci,1)-pts(j,1))/2;
Ay(loci,j+(i-1)*3)=s(ori,3)*qwt(k)*(loc(loci,2)-pts(j,2))/2;
Az(loci,j+(i-1)*3)=s(ori,3)*qwt(k)*(loc(loci,3)-pts(j,3))/2;
Apx(loci,j+(i-1)*3)=s(ori,3)*(loc(loci,1)-pts(j,1))/(2*tvol(i));
Apy(loci,j+(i-1)*3)=s(ori,3)*(loc(loci,2)-pts(j,2))/(2*tvol(i));
Apz(loci,j+(i-1)*3)=s(ori,3)*(loc(loci,3)-pts(j,3))/(2*tvol(i));

    end
    end
end
end
rangec=1:nquad*nt;

%% construct inducance matrix of half-RWG basis functions W=x'*L*x/2 
Lmat=zeros([3*nt 3*nt]);
for j=1:nt
for i=1:nt
Lmat(3*(i-1)+1:3*i,3*(j-1)+1:3*j)=Lmat(3*(i-1)+1:3*i,3*(j-1)+1:3*j)...
    +integralcalcmod(tvol(j),nhat(j,:),...
    p(t2p(i,:),:),p(t2p(j,:),:),del);
end
end

Lmat=10^-7*Lmat;
%change to a loop basis
Atrans=Arwg*Aloop(:,keep);




if nlayers==1
m1=1:3*nt; %layer 1

if debug>=1
Ax1=Ax(:,m1)*Atrans(m1,:);
Ay1=Ay(:,m1)*Atrans(m1,:);
Az1=Az(:,m1)*Atrans(m1,:);
    ntt=min(nt,4);
    for i=1:ntt
        hold on
quiver3(loc(:,1),loc(:,2),loc(:,3),...
    Ax1(:,i),Ay1(:,i),Az1(:,i));
    end
    axis equal
    pause
close all

end
%%find anergy orthonormal modes layer 1
Lmat1=(Atrans(m1,:)'*Lmat(m1,m1)*Atrans(m1,:));
[u1,s1,v1]=svd(Lmat1);
s1(s1~=0)=1./s1(s1~=0);
Lmat1=sqrt(s1)*(v1'*Lmat1*v1)*sqrt(s1);
Atrans(m1,:)=Atrans(m1,:)*v1*sqrt(s1);
Atrans1=v1*sqrt(s1);
Ax1=Ax(:,m1)*Atrans(m1,:);
Ay1=Ay(:,m1)*Atrans(m1,:);
Az1=Az(:,m1)*Atrans(m1,:);

if debug>=1
Ax1=Ax(:,m1)*Atrans(m1,:);
Ay1=Ay(:,m1)*Atrans(m1,:);
Az1=Az(:,m1)*Atrans(m1,:);
    ntt=min(nt,10);
    for i=1:1
        hold on
quiver3(loc(:,1),loc(:,2),loc(:,3),...
    Ax1(:,i),Ay1(:,i),Az1(:,i));
    end
    axis equal
    pause
close all

end
lmax=30;
globind=@(m,l) l+1+m*(lmax-(m-1)/2);

coeff=computemultipolesminenergcoil(loc,2*pi*3000,lmax,.07);
coeff1=coeff(:,:,1)*Ax1+coeff(:,:,2)*Ay1+coeff(:,:,3)*Az1;
[u1 s1 v1]=svd(coeff1);
Atrans1=Atrans1*v1;
nke=24;
Atransp=zeros([numel(keep) nke]);
Atransp(:,1:nke)=Atrans1(:,1:nke);


elseif nlayers==2
m1=1:3*nt/2; %layer 1
Lmat1=(Atrans(m1,:)'*Lmat(m1,m1)*Atrans(m1,:));
%%find energy orthonormal modes layer 1
[u1,s1,v1]=svd(Lmat1);
s1(s1~=0)=1./s1(s1~=0);
Lmat1=sqrt(s1)*(v1'*Lmat1*v1)*sqrt(s1);
Atrans(m1,:)=Atrans(m1,:)*v1*sqrt(s1);

%%find energy orthonormal modes layer 2
m2=3*nt/2+1:3*nt;% layer 2
Lmat2=(Atrans(m2,:)'*Lmat(m2,m2)*Atrans(m2,:));
[u2 s2 v2]=svd(Lmat2);
s2(s2~=0)=1./s2(s2~=0);
Lmat2=sqrt(s2)*(v2'*Lmat2*v2)*sqrt(s2);
Atrans(m2,:)=Atrans(m2,:)*v2*sqrt(s2);
Atrans1=v1*sqrt(s1);
Atrans2=v2*sqrt(s2);
%layer one loops orthonormal
Ax1=Ax(:,m1)*Atrans(m1,:);
Ay1=Ay(:,m1)*Atrans(m1,:);
Az1=Az(:,m1)*Atrans(m1,:);
%layer two loops
Ax2=Ax(:,m2)*Atrans(m2,:);
Ay2=Ay(:,m2)*Atrans(m2,:);
Az2=Az(:,m2)*Atrans(m2,:);
nbas=numel(keep);
lmax=30;
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
coeff=computemultipolesminenergcoil(loc,2*pi*3000,lmax,.07);
coeff1=coeff(:,:,1)*Ax1+coeff(:,:,2)*Ay1+coeff(:,:,3)*Az1;
coeff2=coeff(:,:,1)*Ax2+coeff(:,:,2)*Ay2+coeff(:,:,3)*Az2;
[u1 s1 v1]=svd(coeff1);
[u2 s2 v2]=svd(coeff2);
Atrans1=Atrans1*v1;
Atrans2=Atrans2*v2;

nke=48;
Atransp=zeros([numel(keep) 2*nke]);
Atransp(:,1:2:2*nke)=Atrans1(:,1:nke);
Atransp(:,2:2:2*nke)=Atrans2(:,1:nke);
nke=2*nke;
end
%%
%%
Atrans=Arwg*Aloop(:,keep);
Ax=Ax*Atrans*Atransp;
Ay=Ay*Atrans*Atransp;
Az=Az*Atrans*Atransp;

Apx=Apx*Atrans*Atransp;
Apy=Apy*Atrans*Atransp;
Apz=Apz*Atrans*Atransp;
Lmat2=(Atrans*Atransp)'*Lmat*Atrans*Atransp;
Atrans=Atransp;



if debug>=1
    ntt=min(nt,10);
    for i=1:1
        hold on
quiver3(loc(:,1),loc(:,2),loc(:,3),...
    Ax(:,i),Ay(:,i),Az(:,i));
    end
    axis equal
    pause
close all

end

 %%
 for i=1:nke
     js(:,3)=Az(:,i);
     js(:,2)=Ay(:,i);
     js(:,1)=Ax(:,i);
 harmcoeff{i}=computemultipolesoptimization(loc,js,2*pi*3000,lmax,0.085);
 end
 bigR=.085;
 Asym=Asym(:,keep);

 save(outputfile,'Lmat2','harmcoeff','bigR','lmax','Atrans',...
     'p','t2p','Asym','psym','t2psym','Ax','Ay','Az','Apx','Apy',...
     'Apz','loc','rangec','that1','that2','nhat','-v7.3');
 toc
 %%
 Fsym=scatteredInterpolant(psym(:,1),psym(:,2),psym(:,3),'natural','none');

if debug==1
    i=1;
Att=Asym*Atrans(:,1);
trisurf(t2psym,psym(:,1),psym(:,2),psym(:,3),Att);

ms=max(Att);
mi=min(Att);
del=(ms-mi)/8;
levs=mi+del/2:del:ms-del/2;
[cou{i},hout] = tricontour(psym(:,1:2),t2psym,Att,levs);

close all
[X Y Z]=sphere(20);
X=.085*X;
Y=.085*Y;
Z=.085*Z;

for i=1
    en=0;
    ct=1;
c=cou{i};
while en<numel(c(1,:))
    st=en+1;
en=c(2,st)+st;
if sign(c(1,st))==1
    col='r';
    sig(ct)=1;
else
    sig(ct)=-1;
    col='b';
end
st=st+1;

clear temp

xcoil{ct}=c(1,st:en);
ycoil{ct}=c(2,st:en);
zcoil{ct}=Fsym(c(1,st:en),c(2,st:en));
ct=ct+1;
plot3(c(1,st:en),c(2,st:en),Fsym(c(1,st:en),c(2,st:en)),col);
hold on

end
surf(X,Y,Z,'edgealpha',0)
view([0 0 1]);

end
end