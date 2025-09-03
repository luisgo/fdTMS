harmfile='../harmfillamentary2.mat';
Ncc=22;
%for dept=[49 65 79]
ct=1;
for dept=[100 140 170]

filename='paperover2optimsfilamentary';
filenamesav=strcat(filename,'d',num2str(dept),'ncoils',num2str(Ncc));
designfile=strcat(filenamesav,'lmax30.mat');
ener=[50 100 150 200 250 300];
peakSR2=10^9*ones(size(ener));
%Jmax=10^6*[1.5832    1.9988    3.2843    3.2195    4.0595]*.9;%125
%Jmax=10^6*[3.2052    5.1185    5.6964    7.1927    9.3359]*.9;%96
Jmax=10^12*ones(size(ener));%142
gap=10^30*ones(size(ener));
 %ncoils = runoptiminogaL12layer2(dept,ener(:),Jmax(:),gap(:),Ncc,harmfile,filename);
%ncoils = runoptiminogatwo([dept,dept2(ct)],ener,Jmax,peakSR,Ncc,harmfile,filename);
 ncoils = runoptiminoga2(dept,ener,Jmax,peakSR2,Ncc,harmfile,filename);
 %ncoils = runoptiminoga(dept,ener,Ncc,harmfile,filename);
close all
addpath('C:\Users\ljg24\Desktop\postproc\determinetradeoff\')

addpath('..\determinetradeoff');

[Wi,Si,peakJ,peakSR,Escalp,Escalp12,Ndep,Evol,Edecay,Escalp,Ssq,Ndepsq]=anycoildeterminetradeoff_f(harmfile,designfile,dept,Ncc);
save(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\PAPERFillamentaryWOUT1002',num2str(dept),'.mat'),'Wi','Si','Ssq','peakJ','peakSR','Escalp','Escalp12','Ndep','Evol','Edecay','Escalp','Ndepsq');
peakSR
ct=ct+1;
end
%%

[xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,4,10,Ncc);
Att1=Att1/max(Att1(:));
%%

%%
addpath('E:/Coiloptcode')

Att1=Att1/max(Att1(:));

[xcoilp,ycoilp,zcoilp,zzcoilp]=getcontours(psym,t2psym,Att1,Fsym1,[.5 1.5]/5);
%%
[t2p,pmesh]=TEMPgnerateholes(psym,t2psym,xcoilp,ycoilp,zcoilp,[3 1]);
trisurf(t2p,pmesh(:,1),pmesh(:,2),pmesh(:,3))

hold on

for i=1:numel(xcoil)
    hold on
    if sign(zzcoil{i})==1
plot3(xcoil{i},ycoil{i},zcoil{i},'blue','linewidth',2);
    else
plot3(xcoil{i},ycoil{i},zcoil{i},'red','linewidth',2);
    end
end
%%

%%
for i=1:numel(xcoil)
    len=sqrt((xcoil{i}(2:1:end)-xcoil{i}(1:end-1)).^2+...
             (ycoil{i}(2:1:end)-ycoil{i}(1:end-1)).^2+...
             (zcoil{i}(2:1:end)-zcoil{i}(1:end-1)).^2);


         
         len=cumsum([0;len(:)]);
         len(end)

         N=ceil(len(end)/.0001);
     
    locs=0:len(end)/N:len(end);
     xcoil{i}=interp1(len,xcoil{i},locs,'linear');
     ycoil{i}=interp1(len,ycoil{i},locs,'linear');
     zcoil{i}=interp1(len,zcoil{i},locs,'linear');    
end
save('..\newcoilsp.mat','xcoil','ycoil','zcoil','zzcoil','Fsym1')
i=1;
clear rs js
rs(:,1)=cat(1,(xcoil{i}(2:end)+xcoil{i}(1:end-1))'/2,...
    (xcoil{i}(1)+xcoil{i}(end))/2);
rs(:,2)=cat(1,(ycoil{i}(2:end)+ycoil{i}(1:end-1))'/2,...
    (ycoil{i}(1)+ycoil{i}(end))/2);
rs(:,3)=cat(1,(zcoil{i}(2:end)+zcoil{i}(1:end-1))'/2,...
    (zcoil{i}(1)+zcoil{i}(end))/2);
js(:,1)=cat(1,(xcoil{i}(2:end)-xcoil{i}(1:end-1))',...
    (xcoil{i}(1)-xcoil{i}(end)));
js(:,2)=cat(1,(ycoil{i}(2:end)-ycoil{i}(1:end-1))',...
    (ycoil{i}(1)-ycoil{i}(end)));
js(:,3)=cat(1,(zcoil{i}(2:end)-zcoil{i}(1:end-1))',...
    (zcoil{i}(1)-zcoil{i}(end)));

for i=2:numel(xcoil)
rst(:,1)=cat(1,rs(:,1),(xcoil{i}(2:end)+xcoil{i}(1:end-1))'/2,...
    (xcoil{i}(1)+xcoil{i}(end))/2);
rst(:,2)=cat(1,rs(:,2),(ycoil{i}(2:end)+ycoil{i}(1:end-1))'/2,...
    (ycoil{i}(1)+ycoil{i}(end))/2);
rst(:,3)=cat(1,rs(:,3),(zcoil{i}(2:end)+zcoil{i}(1:end-1))'/2,...
    (zcoil{i}(1)+zcoil{i}(end))/2);
jst(:,1)=cat(1,js(:,1),(xcoil{i}(2:end)-xcoil{i}(1:end-1))',...
    (xcoil{i}(1)-xcoil{i}(end)));
jst(:,2)=cat(1,js(:,2),(ycoil{i}(2:end)-ycoil{i}(1:end-1))',...
    (ycoil{i}(1)-ycoil{i}(end)));
jst(:,3)=cat(1,js(:,3),(zcoil{i}(2:end)-zcoil{i}(1:end-1))',...
    (zcoil{i}(1)-zcoil{i}(end)));
rs=rst;js=jst;
clear rst jst
end
size(rs)
size(js)
addpath('C:\workbackup\ljg24\loopstar\Arraycode\Simtools\InductanceEnergy')
 L=computeinductance(rs,js);
 

quiver3(rs(:,1),rs(:,2),rs(:,3),js(:,1),js(:,2),js(:,3))
%%
%  clear rs js
%   rs=loc;
%        js(:,3)=Az(:,1:Ncc)*Jcalc2{1};
%        js(:,2)=Ay(:,1:Ncc)*Jcalc2{1};
%        js(:,1)=Ax(:,1:Ncc)*Jcalc2{1};
% rv=rs;
% Jdc=js;
addpath('C:\Users\ljg24\Desktop\postproc\lfMFIE\')

[Xeval,Yeval,Zeval]=ndgrid(-.06975:.0005:0,-.06975:.0005:0,.06975:-.0005:0.0345);
n=size(Xeval);
pt(:,1)=Xeval(:);
pt(:,2)=Yeval(:);
pt(:,3)=Zeval(:);
    Efield=computefields(rv,Jdc,zeros(size(rv)),zeros(size(rv(:,1)))...
       ,2*pi*3000,20,.085,pt);
   [Xeval,Yeval,Zeval]=ndgrid(0,0,.07:-.0001:.01);
pt2(:,1)=Xeval(:);
pt2(:,2)=Yeval(:);
pt2(:,3)=Zeval(:);
    Efield2=computefields(rv,Jdc,zeros(size(rv)),zeros(size(rv(:,1)))...
       ,2*pi*3000,20,.085,pt2);
   
rc=sqrt(sum(pt.^2,2));
Efield(rc>.07,:)=0;
plot(Efield2(:,1))

Ema2=sqrt(sum(abs(Efield2).^2,2));
con=50/Ema2(141);
L=computeinductance(rv,Jdc)
ener=L*con^2/2
plot(Ema2)

Ema=con*sqrt(sum(abs(Efield).^2,2));
Ema2=con*sqrt(sum(abs(Efield2).^2,2));



V=sum(Ema(:)>=50)/2;
d=(nnz(Ema2>=50)/10-.1)
S=V*10^-2/d
Ema=reshape(Ema,n);
imagesc(Ema(:,:,round(d)))
%%
load C:\Users\ljg24\Desktop\postproc\generatecoil\firstordercorr.mat nx ny nz

wireth=.00272+.0005;
wirehe=2*.00467;


for i=1:numel(xcoil);
clear rs js
rs(:,1)=cat(1,(xcoil{i}(2:end)+xcoil{i}(1:end-1))'/2,...
    (xcoil{i}(1)+xcoil{i}(end))/2);
rs(:,2)=cat(1,(ycoil{i}(2:end)+ycoil{i}(1:end-1))'/2,...
    (ycoil{i}(1)+ycoil{i}(end))/2);
rs(:,3)=cat(1,(zcoil{i}(2:end)+zcoil{i}(1:end-1))'/2,...
    (zcoil{i}(1)+zcoil{i}(end))/2);
drs=rs(2:end,:)-rs(1:end-1,:);
rs1=(rs(2:end,:)+rs(1:end-1,:))/2;
rs2=(rs(2:end,:)+rs(1:end-1,:))/2;
nh=1./sqrt(sum(drs.^2,2));
drs(:,1)=drs(:,1).*nh;
drs(:,2)=drs(:,2).*nh;
drs(:,3)=drs(:,3).*nh;
clear nhat
nhat(:,1)=nx((rs(2:end,1)+rs(1:end-1,1))/2,(rs(2:end,2)+rs(1:end-1,2))/2);
nhat(:,2)=ny((rs(2:end,1)+rs(1:end-1,1))/2,(rs(2:end,2)+rs(1:end-1,2))/2);
nhat(:,3)=nz((rs(2:end,1)+rs(1:end-1,1))/2,(rs(2:end,2)+rs(1:end-1,2))/2);
tan=cross(drs,nhat);
nh=1./sqrt(sum(tan.^2,2));
tan(:,1)=tan(:,1).*nh;
tan(:,2)=tan(:,2).*nh;
tan(:,3)=tan(:,3).*nh;
rs1=rs1-wireth/2*tan-wirehe*nhat/2;
rs2=rs2+wireth/2*tan-wirehe*nhat/2;
rs3=rs1;rs4=rs2;
rs3=rs3+wirehe*nhat;
rs4=rs4+wirehe*nhat;
addpath('C:\Users\ljg24\Desktop\postproc\generatecoil');
[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);
nte=numel(te2p)/4;
np=numel(p)/3;
if i==1
ntec=nte;npc=np;
te2pc=te2p;pc=p;
portc=port;portidc=portid;
regc=reg;
else
te2pc=cat(1,te2pc,te2p+npc);pc=cat(1,pc,p);
portc=cat(1,portc,port+npc);portidc=cat(1,portidc,portid);
regc=cat(1,regc,reg);
ntec=ntec+nte;npc=npc+np;
end

end
tri=surftri(pc,te2pc);
trisurf(tri,pc(:,1),pc(:,2),pc(:,3))
axis equal
p=pc;te2p=te2pc;reg=regc;port=portc;portid=portidc;
%%
save coilmesh.mat te2p p reg port portid;

addpath('C:\Users\ljg24\Desktop\postproc\lfMFIE\wirecode')

  meshfile='coilmesh.mat'
  [mesh.p,mesh.nte,mesh.nf,mesh.te2f,mesh.te2p,mesh.f2te,mesh.f2p...
 ,mesh.tevol,mesh.tecen,mesh.fvol,mesh.fcen,mesh.reg,mesh.port,mesh.portid,mesh.portf]...
     =createmesh2(meshfile);
 mesh.portcur=zeros(size(mesh.portid));
 x=unique(mesh.portid);
 Nturns=10;
mesh.portid=mesh.portid<=numel(x)/2;
mesh.portid=mesh.portid+1;
tri=surftri(mesh.p,mesh.te2p);
 trisurf(tri,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'edgealpha',0,'facecolor',[127;127;127]/255)
 hold on
  trisurf(mesh.f2p(mesh.portf,:),mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),double(mesh.portid),'edgealpha',0)
  axis equal
  caxis([1 2])
  light

 x=unique(mesh.portid);

 hold on
 axis equal
 

sigma=1;
mesh.portcur(mesh.portid==x(1))=mesh.fvol(mesh.portf(mesh.portid==x(1)))/sum(mesh.fvol(mesh.portf(mesh.portid==x(1))));
mesh.portcur(mesh.portid==x(2))=-mesh.fvol(mesh.portf(mesh.portid==x(2)))/sum(mesh.fvol(mesh.portf(mesh.portid==x(2))));
[Eout,rv,solution]=wirefem_cl(mesh,1,1,sigma);
rv=rv';Eout=Eout';
Eout(:,1)=sigma.*Eout(:,1).*mesh.tevol(:);
Eout(:,2)=sigma.*Eout(:,2).*mesh.tevol(:);
Eout(:,3)=sigma.*Eout(:,3).*mesh.tevol(:);
np=numel(rv(:,1));
%jv(end/2+1:end,:)=-jv(end/2+1:end,:);
quiver3(rv(:,1),rv(:,2),rv(:,3),Eout(:,1),Eout(:,2),Eout(:,3))
jv=2*Nturns*Eout;

vol=mesh.tevol;
L=computeinductance(rv,jv)
p=mesh.p;

[Jdc,rv,Ldc,R]=wirevieLaccurate_cl(meshfile,2*pi*3000,sigma,2*Nturns);

Jdc(:,1)=Jdc(:,1).*mesh.tevol(:);
Jdc(:,2)=Jdc(:,2).*mesh.tevol(:);
Jdc(:,3)=Jdc(:,3).*mesh.tevol(:);