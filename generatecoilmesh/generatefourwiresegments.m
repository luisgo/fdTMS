clear addnum
p1=[0.04 0.01129 .1000];
p2=[0.007672 0.015868 .105];


load ../firstordercorr_middle.mat nx ny nz
h1=-.004067;
Npts=100;
 convi=(1/Npts:1/Npts:1);
 Npts=Npts;
 ax=p1(1)*convi+p2(1)*(1-convi);ay=p1(2)*convi+p2(2)*(1-convi);
 az=p1(3)*convi+p2(3)*(1-convi);
  load ../modeheight.mat
[TH,PHI,R] = cart2sph(-p(:,1),p(:,2),p(:,3));
rad=scatteredInterpolant(TH,PHI,R);
[TH,PHI] = cart2sph(ax,ay,az);
[ax,ay,az]=sph2cart(TH,PHI,rad(TH,PHI));%project
 addnum(1:Npts,1)=ax+h1*nx(ax,ay);
 addnum(1:Npts,2)=ay+h1*ny(ax,ay);
addnum(1:Npts,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum(1:Npts,3)=az+h1*nz(ax,ay);
 addnum(1:Npts,4)=ax;
 addnum(1:Npts,5)=ay;
 wireth=.00575*ones(size(addnum(:,1)));
[rst1,rst2,rst3,rst4,rst]=converttowire(addnum,wireth,wirehe,nx,ny,nz,addnum(:,4),addnum(:,5))

 plot3(rst1(:,1),rst1(:,2),rst1(:,3),'b')
hold on
 plot3(rst2(:,1),rst2(:,2),rst2(:,3),'b')
hold on
 plot3(rst3(:,1),rst3(:,2),rst3(:,3),'b')
hold on
 plot3(rst4(:,1),rst4(:,2),rst4(:,3),'b')
hold on
 plot3(rst(:,1),rst(:,2),rst(:,3),'r')

 [te2pS,pS]=maketetramesh(rst1,rst2,rst3,rst4);
rst1(:,1)=-rst1(:,1);
rst2(:,1)=-rst2(:,1);
rst3(:,1)=-rst3(:,1);
rst4(:,1)=-rst4(:,1);
[te2p1,p1]=maketetramesh(rst1,rst2,rst3,rst4);
 te2pS=cat(1,te2pS,te2p1+numel(pS(:,1)));
 pS=cat(1,pS,p1);
rst1(:,2)=-rst1(:,2);
rst2(:,2)=-rst2(:,2);
rst3(:,2)=-rst3(:,2);
rst4(:,2)=-rst4(:,2);
[te2p1,p1]=maketetramesh(rst1,rst2,rst3,rst4);
 te2pS=cat(1,te2pS,te2p1+numel(pS(:,1)));
 pS=cat(1,pS,p1);
 rst1(:,1)=-rst1(:,1);
rst2(:,1)=-rst2(:,1);
rst3(:,1)=-rst3(:,1);
rst4(:,1)=-rst4(:,1);
[te2p1,p1]=maketetramesh(rst1,rst2,rst3,rst4);
 te2pS=cat(1,te2pS,te2p1+numel(pS(:,1)));
 pS=cat(1,pS,p1);
 
 p1=[0.00413 0.024 .1040];
p2=[0.02 0.03731 .099];
clear addnum
  load ../modeheight.mat

load ../firstordercorr_middle.mat nx ny nz

 ax=p1(1)*convi+p2(1)*(1-convi);ay=p1(2)*convi+p2(2)*(1-convi);
 az=p1(3)*convi+p2(3)*(1-convi);
%   load ../modeheight.mat
% [TH,PHI,R] = cart2sph(-p(:,1),p(:,2),p(:,3));
% rad=scatteredInterpolant(TH,PHI,R);
% [TH,PHI] = cart2sph(ax,ay,az);
% [ax,ay,az]=sph2cart(TH,PHI,rad(TH,PHI));%project
 addnum(1:Npts,1)=ax+h1*nx(ax,ay);
 addnum(1:Npts,2)=ay+h1*ny(ax,ay);
 addnum(1:Npts,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum(1:Npts,4)=ax;
 addnum(1:Npts,5)=ay;
 wireth=.00575*ones(size(addnum(:,1)));
[rst1,rst2,rst3,rst4,rst]=converttowire(addnum,wireth,wirehe,nx,ny,nz,addnum(:,4),addnum(:,5))

 plot3(rst1(:,1),rst1(:,2),rst1(:,3),'b')
hold on
 plot3(rst2(:,1),rst2(:,2),rst2(:,3),'b')
hold on
 plot3(rst3(:,1),rst3(:,2),rst3(:,3),'b')
hold on
 plot3(rst4(:,1),rst4(:,2),rst4(:,3),'b')
hold on
 plot3(rst(:,1),rst(:,2),rst(:,3),'r')

 
  
 [te2p1,p1]=maketetramesh(rst1,rst2,rst3,rst4);
  te2pS=cat(1,te2pS,te2p1+numel(pS(:,1)));
 pS=cat(1,pS,p1);
rst1(:,1)=-rst1(:,1);
rst2(:,1)=-rst2(:,1);
rst3(:,1)=-rst3(:,1);
rst4(:,1)=-rst4(:,1);
[te2p1,p1]=maketetramesh(rst1,rst2,rst3,rst4);
 te2pS=cat(1,te2pS,te2p1+numel(pS(:,1)));
 pS=cat(1,pS,p1);
rst1(:,2)=-rst1(:,2);
rst2(:,2)=-rst2(:,2);
rst3(:,2)=-rst3(:,2);
rst4(:,2)=-rst4(:,2);
[te2p1,p1]=maketetramesh(rst1,rst2,rst3,rst4);
 te2pS=cat(1,te2pS,te2p1+numel(pS(:,1)));
 pS=cat(1,pS,p1);
 rst1(:,1)=-rst1(:,1);
rst2(:,1)=-rst2(:,1);
rst3(:,1)=-rst3(:,1);
rst4(:,1)=-rst4(:,1);
[te2p1,p1]=maketetramesh(rst1,rst2,rst3,rst4);
 te2pS=cat(1,te2pS,te2p1+numel(pS(:,1)));
 pS=cat(1,pS,p1);
tri=surftri(pS,te2pS);
 save fourcorners.mat tri pS
 