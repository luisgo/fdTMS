clear all
addpath('C:\Users\ljg24\Desktop\shamcoil\optimize');
load coildes.mat xcoil ycoil zcoil Fsym1;
xcoiltemp=xcoil; ycoiltemp=ycoil; zcoiltemp=zcoil;
xcoil2=xcoil; ycoil2=ycoil; zcoil2=zcoil;
for i=1:numel(xcoil);
rct(i,1)=mean(xcoil{i});
rct(i,2)=mean(ycoil{i});
rct(i,3)=mean(zcoil{i});
nctr(i)=numel(xcoil{i});
end
load ../newcoilsp.mat xcoil ycoil zcoil Fsym1;
for i=1:numel(xcoil);
rct2(i,1)=mean(xcoil{i});
rct2(i,2)=mean(ycoil{i});
rct2(i,3)=mean(zcoil{i});
nctr2(i)=numel(xcoil{i});

end
%%
Amat=zeros([numel(xcoil) 1])
ct=zeros([numel(xcoiltemp) 1]);

    for j=1:numel(xcoiltemp) 
for i=1:numel(xcoil)
        Amat(i)=norm(rct(j,:)-rct2(i,:));
end
if j>3 && j~=5 && j~=7
    [B,I]=sort(Amat);
I=I(B(:)<2*B(1));
[B,I2]=sort(abs(nctr(j)-nctr2(I)));
 B(1)
I=I(I2);
    ct(i)=I(1);

else
    
kk=[1 2 3 0 5 0 8];
    I(1)=kk(j);
    ct(i)=I(1);
end


[B,I]=sort(sqrt((xcoiltemp{j}(1)-xcoil{ct(i)}(:)).^2+...
                (ycoiltemp{j}(1)-ycoil{ct(i)}(:)).^2+...
                     (zcoiltemp{j}(1)-zcoil{ct(i)}(:)).^2));
    rval=sqrt((rct2(ct(i),1)-xcoil{ct(i)}(:)).^2+...
              (rct2(ct(i),2)-ycoil{ct(i)}(:)).^2);
    anglesx=atan2((ycoil{ct(i)}(:)-rct2(ct(i),2)),(xcoil{ct(i)}(:)-rct2(ct(i),1)));

rval=sqrt((rct(j,1)-xcoiltemp{j}(:)).^2+...
          (rct(j,2)-ycoiltemp{j}(:)).^2);

      anglesxtemp=atan2((ycoiltemp{j}(:)-rct(j,2)),(xcoiltemp{j}(:)-rct(j,1)));
    
    [B,I]=sort(abs(anglesxtemp(1)-anglesx));
plot(anglesx)
            if I(1)==1
    xcoiltemp{j}=xcoil{ct(i)};
    ycoiltemp{j}=ycoil{ct(i)};
    zcoiltemp{j}=zcoil{ct(i)};
            else
    xcoiltemp{j}=cat(2,xcoil{ct(i)}(I(1):end),xcoil{ct(i)}(1:I(1)-1));
    ycoiltemp{j}=cat(2,ycoil{ct(i)}(I(1):end),ycoil{ct(i)}(1:I(1)-1));
    zcoiltemp{j}=cat(2,zcoil{ct(i)}(I(1):end),zcoil{ct(i)}(1:I(1)-1));
            end

    end
xcoil=xcoiltemp;
ycoil=ycoiltemp;
zcoil=zcoiltemp;


sqrt((xcoil{2}(2:end)-xcoil{2}(1:end-1)).^2+...
     (ycoil{2}(2:end)-ycoil{2}(1:end-1)).^2+...
     (zcoil{2}(2:end)-zcoil{2}(1:end-1)).^2);
load ../firstordercorr_middle.mat nx ny nz
%%
subplot(1,2,1),
for i=1:numel(xcoil)
plot3(xcoil{i},ycoil{i},1:numel(zcoil{i}))
hold on
text(xcoil{i}(1),ycoil{i}(1),1,num2str(i));
% for j=1:200:numel(xcoil{i})
view(2)
% text(xcoil{i}(j),ycoil{i}(j),zcoil{i}(j),num2str(j));
% end
end
subplot(1,2,2),
for i=1:numel(xcoil2)
plot3(xcoil2{i},ycoil2{i},1:numel(zcoil2{i}))
hold on
text(xcoil2{i}(1),ycoil2{i}(1),1,num2str(i));
% for j=1:200:numel(xcoil{i})
% text(xcoil{i}(j),ycoil{i}(j),zcoil{i}(j),num2str(j));
% end
end
view(2)
%%
%%inner middle coil
h1=-.004067;
h2=0;
c1=[1,2,3,5,7];
%side L
c2=[8];
%side R
c3=[9];
%bottom
c3=[4,6];
res=.0001


layerp=h1;
save layerinfo.mat layerp;

rs=coilsegment(xcoil,ycoil,zcoil,8,[235:448],h1,h1,[],res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);


%rs=coilsegment(xcoil,ycoil,zcoil,8,[235:numel(xcoil{8}) 1:448],h2,h1,rs,res,Fsym1...
%    ,xcoil2,ycoil2,zcoil2);

%%xx



 rs=coilsegment(xcoil,ycoil,zcoil,7,[214:419],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);

rs=coilsegment(xcoil,ycoil,zcoil,5,[378:numel(xcoil{5}) 1:337],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);

rs=coilsegment(xcoil,ycoil,zcoil,3,[261:numel(xcoil{3}) 1:220],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
rs=coilsegment(xcoil,ycoil,zcoil,2,[230:numel(xcoil{2}) 1:189],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
%rs=coilsegment(xcoil,ycoil,zcoil,1,[115:numel(xcoil{1}) 1:74],h1,h1,rs,res,Fsym1...
%    ,xcoil2,ycoil2,zcoil2);
 addnum=rs(end,:);
 convi=(1/30:1/30:1);
 
  ax=rs(end,1)*(1-convi);ay=rs(end,2)-.005*convi;
 addnum(1:30,1)=ax+h1*nx(ax,ay);
 addnum(1:30,2)=ay+h1*ny(ax,ay);
 addnum(1:30,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum(1:30,4)=ax;
 addnum(1:30,5)=ay;
 
 rs=cat(1,rs,addnum);
 rs=coilsegment(xcoil,ycoil,zcoil,1,[115:numel(xcoil{1}) 1:74],h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
 rs=coilsegment(xcoil,ycoil,zcoil,2,[230:numel(xcoil{2}) 1:189],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
 rs=coilsegment(xcoil,ycoil,zcoil,3,[261:numel(xcoil{3}) 1:220],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
 rs=coilsegment(xcoil,ycoil,zcoil,5,[378:numel(xcoil{5}) 1:337],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
%feed 1
rs=coilsegment(xcoil,ycoil,zcoil,7,[436:683],h1,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
rs=coilsegment(xcoil,ycoil,zcoil,9,[221:numel(xcoil{9}) 1:191] ,h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
 rs=coilsegment(xcoil,ycoil,zcoil,9,[215:numel(xcoil{9}) 1:140] ,h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
 rs=coilsegment(xcoil,ycoil,zcoil,7,[763:numel(xcoil{7}) 1:150],h1,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
 rs=coilsegment(xcoil,ycoil,zcoil,7,[763:numel(xcoil{7}) 1:150],h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);


Npts=50;
plot3(rs(:,1),rs(:,2),rs(:,3));
rs=[];

%% 

rs=coilsegment(xcoil,ycoil,zcoil,7,[460:numel(xcoil{7}) 1:439],h1,h2,rs,res,Fsym1...
     ,xcoil2,ycoil2,zcoil2);

 rs=coilsegment(xcoil,ycoil,zcoil,6,[4606:numel(xcoil{6}) 1:4585] ,h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
addnum=rs(end,:);
Npts=50;

rs=coilsegment(xcoil,ycoil,zcoil,6,[4626:numel(xcoil{6}) 1:1473] ,h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);

convi=(1/Npts:1/Npts:1);

 ax=rs(1,1)*(1-convi);ay=rs(1,2)+.005*convi;
 addnum(1:Npts,1)=ax+h1*nx(ax,ay);
 addnum(1:Npts,2)=ay+h1*ny(ax,ay);
 addnum(1:Npts,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum(1:Npts,4)=ax;
 addnum(1:Npts,5)=ay;
 stid=5;
 rs=cat(1,addnum,rs);


 % %round corners
rs=coilsegment(xcoil,ycoil,zcoil,4,[673:numel(xcoil{4}) 1:638] ,h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
rs=coilsegment(xcoil,ycoil,zcoil,4,[628:numel(xcoil{4}) 1:590] ,h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
rs=coilsegment(xcoil,ycoil,zcoil,6,[1564:4626],h1,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2);
%y=29 to 32
 % 
 
%%

 convi=(1/Npts:1/Npts:1);
  ax=rs(end,1)*(1-convi);ay=rs(end,2)+.005*convi;
 addnum(1:Npts,1)=ax+h1*nx(ax,ay);
 addnum(1:Npts,2)=ay+h1*ny(ax,ay);
 addnum(1:Npts,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum(1:Npts,4)=ax;
 addnum(1:Npts,5)=ay;
 stid=numel(rs(:,1));
 rs=cat(1,rs,addnum);
Npad=40;
rs(stid-Npad:stid+Npts,1)= movmean(rs(stid-Npad:stid+Npts,1),Npad);
rs(stid-Npad:stid+Npts,2)= movmean(rs(stid-Npad:stid+Npts,2),Npad);
rs(stid-Npad:stid+Npts,3)= movmean(rs(stid-Npad:stid+Npts,3),Npad);
rs(stid-Npad:stid+Npts,4)= movmean(rs(stid-Npad:stid+Npts,4),Npad);
rs(stid-Npad:stid+Npts,5)= movmean(rs(stid-Npad:stid+Npts,5),Npad);

%%%%%%%%%%extra smoothing
Npad=80;
stid=12000;
Npts=1000;
rs(stid-Npad:stid+Npts,1)= movmean(rs(stid-Npad:stid+Npts,1),Npad);
rs(stid-Npad:stid+Npts,2)= movmean(rs(stid-Npad:stid+Npts,2),Npad);
rs(stid-Npad:stid+Npts,3)= movmean(rs(stid-Npad:stid+Npts,3),Npad);
rs(stid-Npad:stid+Npts,4)= movmean(rs(stid-Npad:stid+Npts,4),Npad);
rs(stid-Npad:stid+Npts,5)= movmean(rs(stid-Npad:stid+Npts,5),Npad);
stid=28500;
Npts=1000;
rs(stid-Npad:stid+Npts,1)= movmean(rs(stid-Npad:stid+Npts,1),Npad);
rs(stid-Npad:stid+Npts,2)= movmean(rs(stid-Npad:stid+Npts,2),Npad);
rs(stid-Npad:stid+Npts,3)= movmean(rs(stid-Npad:stid+Npts,3),Npad);
rs(stid-Npad:stid+Npts,4)= movmean(rs(stid-Npad:stid+Npts,4),Npad);
rs(stid-Npad:stid+Npts,5)= movmean(rs(stid-Npad:stid+Npts,5),Npad);

 %%
 
    len=sqrt((rs(2:end,1)-rs(1:end-1,1)).^2+...
             (rs(2:end,2)-rs(1:end-1,2)).^2+...
             (rs(2:end,3)-rs(1:end-1,3)).^2);
          [B,I]=sort(len);
          I=I(B==0);
          rs(I,:)=[];
              len=sqrt((rs(2:end,1)-rs(1:end-1,1)).^2+...
             (rs(2:end,2)-rs(1:end-1,2)).^2+...
             (rs(2:end,3)-rs(1:end-1,3)).^2);
         len=cumsum([0;len(:)]);
    N=ceil(len(end)/.001);
    locs=0:len(end)/N:len(end);
    
     xcoil2=interp1(len,rs(1:end,1),locs,'line');
     ycoil2=interp1(len,rs(1:end,2),locs,'line');
     zcoil2=interp1(len,rs(1:end,3),locs,'line');
     xh=interp1(len,rs(1:end,4),locs,'line');
     yh=interp1(len,rs(1:end,5),locs,'line');
clear rs
 rs(:,1)=xcoil2;
 rs(:,2)=ycoil2;
 rs(:,3)=zcoil2;
save sendtokop.mat rs;

wireth=.0028%+.0005;
wirehe=.00467;
drs=rs(2:end,:)-rs(1:end-1,:);
rs1=(rs(2:end,:)+rs(1:end-1,:))/2;
rs2=(rs(2:end,:)+rs(1:end-1,:))/2;
nh=1./sqrt(sum(drs.^2,2));
drs(:,1)=drs(:,1).*nh;
drs(:,2)=drs(:,2).*nh;
drs(:,3)=drs(:,3).*nh;
clear nhat
nhat(:,1)=nx((xh(2:end)+xh(1:end-1))/2,(yh(2:end)+yh(1:end-1))/2);
nhat(:,2)=ny((xh(2:end)+xh(1:end-1))/2,(yh(2:end)+yh(1:end-1))/2);
nhat(:,3)=nz((xh(2:end)+xh(1:end-1))/2,(yh(2:end)+yh(1:end-1))/2);
tan=cross(drs,nhat);
nh=1./sqrt(sum(tan.^2,2));
tan(:,1)=tan(:,1).*nh;
tan(:,2)=tan(:,2).*nh;
tan(:,3)=tan(:,3).*nh;
rs=rs1+wirehe*nhat/2;
rs1=rs1-wireth/2*tan;
rs2=rs2+wireth/2*tan;
rs3=rs1;rs4=rs2;
rs3=rs3+wirehe*nhat;
rs4=rs4+wirehe*nhat;
save sendtokop.mat rs;

[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);

%%
 plot3(rs1(:,1),rs1(:,2),rs1(:,3),'b')
hold on
 plot3(rs2(:,1),rs2(:,2),rs2(:,3),'b')
hold on
 plot3(rs3(:,1),rs3(:,2),rs3(:,3),'b')
hold on
 plot3(rs4(:,1),rs4(:,2),rs4(:,3),'b')
hold on
 plot3(rs(:,1),rs(:,2),rs(:,3),'r')
%%
nl=numel(p(:,1));
p2=p;
p2(:,2)=-p2(:,2);
te2p2=te2p+nl;
p=cat(1,p,p2);
te2p=cat(1,te2p,te2p2);
reg=cat(1,reg(:),reg(:));
portid=[1,1,2,2,1,1,2,2];
port(end+1:end+4,:)=port(1:4,:)+nl;
save coilmesh.mat te2p p reg port portid;
t2p=surftri(p,te2p);
  tr.faces=t2p;
  tr.vertices=p;
stlwrite('C:\workbackup\ljg24\loopstar\Arraycode\coilwindingsrena.stl',tr)
addpath('C:\Users\ljg24\Desktop\optimisationcodes\anycoils\lfMFIE\wirecode')
addpath('C:\Users\ljg24\Desktop\optimisationcodes\anycoils\')

%%


%%
% 
% [p, te2p ]=tet_mesh_refine (p,te2p);
% %%
% close all
% t2p=surftri(p,te2p);
% trisurf(t2p,p(:,1),p(:,2),p(:,3))
% pcen=(p(t2p(:,1),:)+p(t2p(:,2),:)+p(t2p(:,3),:))/3;
% for i=1:numel(t2p)
% text(pcen(i,1),pcen(i,2),pcen(i,3),num2str(i));
% end
%%
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
 

sigma=1;
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
jv=2*Nturns*Eout;

vol=mesh.tevol;
L=computeinductance(rv,jv)
p=mesh.p;

%%
[Jdc,rv,Ldc,R]=wirevieLaccurate_cl(meshfile,2*pi*3000,sigma,2*Nturns);
Jdc(:,1)=Jdc(:,1).*mesh.tevol(:);
Jdc(:,2)=Jdc(:,2).*mesh.tevol(:);
Jdc(:,3)=Jdc(:,3).*mesh.tevol(:);
   norm(Jdc(:)+jv(:))/norm(Jdc(:))
   
   
save fdtmscoilsp.mat rv jv tri p Jdc R Ldc;