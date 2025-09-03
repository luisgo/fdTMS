clear all
[47 62 63 71]
harmfile='../newwireharmonics.mat';
dept=71;Wc=107;Sc=15.2;
Ncc=24;
designfile=strcat('optimT3no707d',num2str(dept),'ncoils24lmax30.mat');
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    ct=1;
    for i=6-numel(Wi):5
[xcoilc{ct},ycoilc{ct},zcoilc{ct},zzcoilc{ct},Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,i,10,Ncc);
Att1=Att1/max(Att1(:));
ct=ct+1;
end
close all

for i=1:numel(Wi)
    subplot(2,3,i),
    title(strcat('(S,W)=(',num2str(round(10*Si(i))/10),'cm^2,',num2str(round(Wi(i))),'J)'));
    for j=1:numel(xcoilc{i})
%for j=[1 2 3 5 8 12 15 16 17 18]
%for j=[4 6 10 14]
%for j=[7 9 11 13]
    hold on
    if sign(zzcoilc{i}{j})==1
plot3(xcoilc{i}{j},ycoilc{i}{j},zcoilc{i}{j},'blue','linewidth',2);
    else
plot3(xcoilc{i}{j},ycoilc{i}{j},zcoilc{i}{j},'red','linewidth',2);
    end
   % text(xcoilc{i}{j}(1),ycoilc{i}{j}(1),zcoilc{i}{j}(1),num2str(j));
    
       end
    axis equal
end


%%

for ival=1:numel(xcoilc)
xcoil=xcoilc{ival};
ycoil=ycoilc{ival};
zcoil=zcoilc{ival};
zzcoil=zzcoilc{ival};

clen=0;
for i=1:numel(xcoil)
    len=sqrt((xcoil{i}(2:1:end)-xcoil{i}(1:end-1)).^2+...
             (ycoil{i}(2:1:end)-ycoil{i}(1:end-1)).^2+...
             (zcoil{i}(2:1:end)-zcoil{i}(1:end-1)).^2);


         
         len=cumsum([0;len(:)]);
        clen=clen+len(end);

         N=ceil(len(end)/.0008);
     
    locs=0:len(end)/N:len(end);
     xcoil{i}=interp1(len,xcoil{i},locs,'linear');
     ycoil{i}=interp1(len,ycoil{i},locs,'linear');
     zcoil{i}=interp1(len,zcoil{i},locs,'linear');    
end
%save('..\newcoilsp.mat','xcoil','ycoil','zcoil','zzcoil','Fsym1')
%i=1;
%i=4;
i=7;
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
%for i=[2 3 5 8 12 15 16 17 18]
%for i=[6 10 14]
%for i=[9 11 13]

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
rs2=rs;
rs2(:,3)=rs(:,3)+.002;
rs(:,3)=rs(:,3)-.002;
js2=js;
rs=cat(1,rs,rs2);
js=cat(1,js,js2);
addpath('C:\workbackup\ljg24\loopstar\Arraycode\Simtools\InductanceEnergy')
 Ldc(ival)=computeinductance(rs,js)
 Rdc(ival)=2*clen/(4*2*58.14)
%plansinglewinding('..\newcoilsp.mat','tryout')
end
Rdc
Ldc
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