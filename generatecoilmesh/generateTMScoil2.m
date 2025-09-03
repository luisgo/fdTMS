function =figure8coilmodel(Nloop,ID,OD,wirehe,wireth);

addpath('C:\Users\ljg24\Desktop\optimisationcodes\anycoils\lfMFIE\wirecode')
addpath('C:\Users\ljg24\Desktop\optimisationcodes\anycoils\')
Nloop=15;
ID=.02;
OD=.072;
NN=1000;
wirehe=.006;
wireth=.0015;
rs(:,1)=linspace(ID/2+wireth/2, OD/2-wireth/2, NN*Nloop).*cos(2*pi*(0:(Nloop*NN-1))/NN)-OD/2+wireth/2;
rs(:,2)=linspace(ID/2+wireth/2, OD/2-wireth/2, NN*Nloop).*sin(2*pi*(0:(Nloop*NN-1))/NN);
rs(:,3)=0;
rs(end+1:2*end,1)=-linspace(OD/2-wireth/2, ID/2+wireth/2, NN*Nloop).*cos(2*pi*(0:(Nloop*NN-1))/NN)+OD/2-wireth/2;
rs(NN*Nloop+1:end,2)=linspace(OD/2-wireth/2, ID/2+wireth/2, NN*Nloop).*sin(2*pi*(0:(Nloop*NN-1))/NN);
rs(NN*Nloop+1:end,3)=0;
rs(end-100:end,3)=linspace(0,wirehe+.00019,101);
 m=numel(rs)/3;
 rs(end+1:end+1000,1)=linspace(rs(end,1),rs(1,1),1000);
  rs(m+1:m+1000,2)=rs(m,2);
  rs(m+1:m+1000,3)=wirehe+.0002;

len=sqrt((rs(2:end,1)-rs(1:end-1,1)).^2+...
         (rs(2:end,2)-rs(1:end-1,2)).^2+...
         (rs(2:end,3)-rs(1:end-1,3)).^2);
         len=cumsum([0;len(:)]);
    N=ceil(len(end)/.001);
    locs=0:len(end)/N:len(end);
    
     xcoil2=interp1(len,rs(1:end,1),locs,'line');
     ycoil2=interp1(len,rs(1:end,2),locs,'line');
     zcoil2=interp1(len,rs(1:end,3),locs,'line');

     %
axis equal
clear rs
 rs(:,1)=xcoil2;
 rs(:,2)=ycoil2;
 rs(:,3)=zcoil2;

rs1=(rs(2:end,:)+rs(1:end-1,:))/2;
rs2=(rs(2:end,:)+rs(1:end-1,:))/2;

drs=rs(2:end,:)-rs(1:end-1,:);
%tan=drs(2:end,:)-drs(1:end-1,:);
%tan(end+1,:)=tan(end,:);
nhat=zeros(size(drs));
nhat(:,3)=1;
tan=cross(drs,nhat);
nh=1./sqrt(sum(drs.^2,2));
drs(:,1)=drs(:,1).*nh;
drs(:,2)=drs(:,2).*nh;
drs(:,3)=drs(:,3).*nh;
tan(:,3)=0;
nh=sqrt(sum(tan.^2,2));
nh=sign(nhat(:,3))./sqrt(sum(tan.^2,2));
tan(:,1)=tan(:,1).*nh;
tan(:,2)=tan(:,2).*nh;
tan(:,3)=tan(:,3).*nh;
nh=sign(nhat(:,3))./sqrt(sum(nhat.^2,2));
nhat(:,1)=nhat(:,1).*nh;
nhat(:,2)=nhat(:,2).*nh;
nhat(:,3)=nhat(:,3).*nh;

rs=rs1+wirehe*nhat/2;
ct=0;


rs1(:,1)=rs1(:,1)-wireth.*tan(:,1)/2;
rs1(:,2)=rs1(:,2)-wireth.*tan(:,2)/2;
rs1(:,3)=rs1(:,3)-wireth.*tan(:,3)/2;
rs2(:,1)=rs2(:,1)+wireth.*tan(:,1)/2;
rs2(:,2)=rs2(:,2)+wireth.*tan(:,2)/2;
rs2(:,3)=rs2(:,3)+wireth.*tan(:,3)/2;
rs3=rs1;rs4=rs2;
rs3=rs3+wirehe*nhat;
rs4=rs4+wirehe*nhat;

%plot3(rs1(:,1),rs1(:,2),rs1(:,3),'b')
%hold on
%plot3(rs2(:,1),rs2(:,2),rs2(:,3),'r')
%plot3(rs3(:,1),rs3(:,2),rs3(:,3),'g')
%plot3(rs4(:,1),rs4(:,2),rs4(:,3),'black')

[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);
 save coilmesh.mat te2p p reg port portid;
t2p=surftri(p,te2p);
trisurf(t2p,p(:,1),p(:,2),p(:,3));

axis equal
  %%%%%%generate coil model
  meshfile='coilmesh.mat'
  [mesh.p,mesh.nte,mesh.nf,mesh.te2f,mesh.te2p,mesh.f2te,mesh.f2p...
 ,mesh.tevol,mesh.tecen,mesh.fvol,mesh.fcen,mesh.reg,mesh.port,mesh.portid,mesh.portf]...
     =createmesh2(meshfile);
 mesh.portcur=zeros(size(mesh.portid));
 x=unique(mesh.portid);
 Nturns=1
mesh.portid=mesh.portid<=numel(x)/2;
mesh.portid=mesh.portid+1;
tri=surftri(mesh.p,mesh.te2p);
 trisurf(tri,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'edgealpha',0,'facecolor',[127;127;127]/255)
 hold on
  trisurf(mesh.f2p(mesh.portf,:),mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),double(mesh.portid),'edgealpha',0)
  axis equal
  caxis([1 2])
  light
%%
 x=unique(mesh.portid);

 hold on
 axis equal
 

sigma=5.96*10^7;
mesh.portcur(mesh.portid==x(1))=mesh.fvol(mesh.portf(mesh.portid==x(1)))/sum(mesh.fvol(mesh.portf(mesh.portid==x(1))));
mesh.portcur(mesh.portid==x(2))=-mesh.fvol(mesh.portf(mesh.portid==x(2)))/sum(mesh.fvol(mesh.portf(mesh.portid==x(2))));
[Eout,rv,solution]=wirefem_cl(mesh,1,1,sigma);
rv=rv';Eout=Eout';
Eout(:,1)=sigma*Eout(:,1).*mesh.tevol(:);
Eout(:,2)=sigma*Eout(:,2).*mesh.tevol(:);
Eout(:,3)=sigma*Eout(:,3).*mesh.tevol(:);
np=numel(rv(:,1));
%jv(end/2+1:end,:)=-jv(end/2+1:end,:);
quiver3(rv(:,1),rv(:,2),rv(:,3),Eout(:,1),Eout(:,2),Eout(:,3))
jv=Eout;

vol=mesh.tevol;
L=computeinductance(rv,jv)
p=mesh.p;

%%
[Jdc,rv,Ldc,R]=wirevieLaccurate_cl(meshfile,2*pi*3000,sigma,Nturns);
Jdc(:,1)=Jdc(:,1).*mesh.tevol(:);
Jdc(:,2)=Jdc(:,2).*mesh.tevol(:);
Jdc(:,3)=Jdc(:,3).*mesh.tevol(:);
   norm(Jdc(:)+jv(:))/norm(Jdc(:))
   
   
