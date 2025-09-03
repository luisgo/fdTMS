function separateloops(xcoilbig,ycoilbig,zcoilbig,...
          xcoilfig8,ycoilfig8,zcoilfig8,...
          xcoilside,ycoilside,zcoilside,outfile);
      res=.0001;
load ../generatecoilmesh/coildes.mat xcoil ycoil zcoil Fsym1;
xcoiltemp=xcoil; ycoiltemp=ycoil; zcoiltemp=zcoil;
xcoil2=xcoil; ycoil2=ycoil; zcoil2=zcoil;
for i=1:numel(xcoil);
rct(i,1)=mean(xcoil{i});
rct(i,2)=mean(ycoil{i});
rct(i,3)=mean(zcoil{i});
nctr(i)=numel(xcoil{i});
end
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
[xcoilbig,ycoilbig,zcoilbig]=removenegside(xcoilbig,ycoilbig,zcoilbig);
[xcoilfig8,ycoilfig8,zcoilfig8]=removenegside(xcoilfig8,ycoilfig8,zcoilfig8);
[xcoilside,ycoilside,zcoilside]=removenegside(xcoilside,ycoilside,zcoilside);
[xcoilbig,ycoilbig,zcoilbig]=reorientcoil(xcoilbig,ycoilbig,zcoilbig);
[xcoilfig8,ycoilfig8,zcoilfig8]=reorientcoil(xcoilfig8,ycoilfig8,zcoilfig8);
[xcoilside,ycoilside,zcoilside]=reorientcoil(xcoilside,ycoilside,zcoilside);
load firstordercorr_middle.mat nx ny nz
for i=1:numel(zcoilbig)
   [~,begbig(i)]=min(1000*(xcoilbig{i}-.008).^2+...
                     (ycoilbig{i}-max(ycoilbig{i})).^2+...
                     (zcoilbig{i}-min(zcoilbig{i})).^2);
  
   [~,endbig(i)]=min(1000*(xcoilbig{i}-.001).^2+...
                     (ycoilbig{i}-max(ycoilbig{i})).^2+...
                     (zcoilbig{i}-min(zcoilbig{i})).^2);
end

 wid1=[];
 wid2=[];
layerp=h1;
save(strcat(outfile,'_layerinfo.mat'),'layerp');
nrs=100000;
[~,stval]=min(xcoilbig{numel(zcoilbig)});
stval=stval-400
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilbig,ycoilbig,zcoilbig,numel(zcoilbig),[stval:numel(xcoil{numel(zcoilbig)}) 1:endbig(numel(zcoilbig))] ,h1,h1,[],res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,6,nrs);
%feed 1

% 
% Npts=50;
% wid2(1)=wid2(1)+50;
%  convi=(1/Npts:1/Npts:1-1/Npts);
%  Npts=Npts-1;
%  ax=rs(1,1)*ones(size(convi));ay=rs(1,2)*(convi);
%  az=rs(1,3)*ones(size(convi));
%   load ../modeheight.mat
% [TH,PHI,R] = cart2sph(p(:,1),-p(:,2),p(:,3));
% rad=scatteredInterpolant(TH,PHI,R);
% [TH,PHI] = cart2sph(ax,ay,az);
% [ax,ay,az]=sph2cart(TH,PHI,rad(TH,PHI));%project
%  addnum(1:Npts,1)=ax+h1*nx(ax,ay);
%  addnum(1:Npts,2)=ay+h1*ny(ax,ay);
%  %addnum(1:Npts,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
%  addnum(1:Npts,3)=az+h1*nz(ax,ay);
%  addnum(1:Npts,4)=ax;
%  addnum(1:Npts,5)=ay;
%  stid=5;
% 
%   addnum2(1:Npts,1)=ax+h2*nx(ax,ay);
%  addnum2(1:Npts,2)=ay+h2*ny(ax,ay);
%  %addnum(1:Npts,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
%  addnum2(1:Npts,3)=az+h2*nz(ax,ay);
%  addnum2(1:Npts,4)=ax;
%  addnum2(1:Npts,5)=ay;
%  addnum2=addnum2(end:-1:1,:);
% 
%  plot3(addnum(:,1),addnum(:,2),addnum(:,3))
%  hold on
%  
%  plot3(addnum2(:,1),addnum2(:,2),addnum2(:,3),'--r')
%  axis equal
%  view(2)
%  rs=cat(1,addnum(1:end-5,:),rs(5:end,:));
%%bigloops
 [rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment2(xcoilbig,ycoilbig,zcoilbig,numel(zcoilbig)-1,[begbig(numel(zcoilbig)-1):numel(xcoilbig{numel(zcoilbig)-1}) 1:endbig(numel(zcoilbig)-1)] ,h1,h1,rs,res,Fsym1,nrs);
%[rs,nrs]=coilsegment2(xcoilbig,ycoilbig,zcoilbig,numel(zcoilbig)-2,[begbig(numel(zcoilbig)-2):numel(xcoilbig{numel(zcoilbig)-2}) 1:endbig(numel(zcoilbig)-2)] ,h1,h1,rs,res,Fsym1...
%    ,nrs);
%[rs,nrs]=coilsegment2(xcoilbig,ycoilbig,zcoilbig,numel(zcoilbig)-3,[451:numel(xcoilbig{numel(zcoilbig)-3}) 1:400] ,h1,h1,rs,res,Fsym1...
%    ,nrs);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilbig,ycoilbig,zcoilbig,numel(zcoilbig),[begbig(numel(zcoilbig)):4566],h1,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,6,nrs);
close all
plot3(rs(:,1),rs(:,2),rs(:,3));
hold on
for i=1:numel(xcoilbig)
plot3(xcoilbig{i},ycoilbig{i},1:numel(xcoilbig{i}));
end

%%fig81andside
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,5,[455:703],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,7,nrs);
nnn=numel(rs(:,1));
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilside,ycoilside,zcoilside,2,[221:numel(xcoilside{2}) 1:171] ,h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,9,nrs);

 [rs,nrs]=coilsegment(xcoilside,ycoilside,zcoilside,4,[215:numel(xcoilside{4}) 1:120] ,h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,9,nrs);
nrstore=nrs;
rsside=rs(nnn+1:end,:);
rsside(:,1)=-rsside(:,1);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,5,[783:numel(xcoilfig8{5}) 1:100],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,7,nrs);
nrs=nrstore;
rs=cat(1,rs,rsside(end:-1:1,:));


%%%sideloop1 
 %[~,~,wid1(end+1),wid2(end+1)]=coilsegment(xcoilside,ycoilside,zcoilside,1,[579:numel(xcoilside{1}) 1:509],h1,h1,rs,res,Fsym1...
 %   ,xcoil2,ycoil2,zcoil2,8,nrs);
%[rs,nrs]=coilsegment(xcoilside,ycoilside,zcoilside,3,[535:numel(xcoilside{3}) 1:448],h1,h1,rs,res,Fsym1...
%   ,xcoil2,ycoil2,zcoil2,8,nrs);

 [rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,5,[214:419],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,7,nrs);

[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,4,[378:numel(xcoilfig8{4}) 1:337],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,5,nrs);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,3,[261:numel(xcoilfig8{3}) 1:220],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,3,nrs);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,2,[230:numel(xcoilfig8{2}) 1:189],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,2,nrs);

 addnum=rs(end,:);
 convi=(1/30:1/30:1);

  ax=rs(end,1)*(1-convi);ay=rs(end,2)-.005*convi;
 addnum(1:30,1)=ax+h1*nx(ax,ay);
 addnum(1:30,2)=ay+h1*ny(ax,ay);
 addnum(1:30,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum(1:30,4)=ax;
 addnum(1:30,5)=ay;
 
 rs=cat(1,rs,addnum);
 
 [rs,nrs]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,1,[115:numel(xcoilfig8{1}) 1:95],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,1,nrs);

 [rs,nrs]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,1,[115:numel(xcoilfig8{1}) 1:74],h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,1,nrs);
 [rs,nrs]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,2,[230:numel(xcoilfig8{2}) 1:189],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,2,nrs);
 [rs,nrs]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,3,[261:numel(xcoilfig8{3}) 1:220],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,3,nrs);
 [rs,nrs]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,4,[378:numel(xcoilfig8{4}) 1:337],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,5,nrs);
 [rs,nrs]=coilsegment(xcoilfig8,ycoilfig8,zcoilfig8,5,[460:numel(xcoilfig8{5}) 1:429],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,7,nrs);
[~,stind]=min((xcoilbig{2}-rs(end,1)).^2+...
              (ycoilbig{2}-rs(end,2)).^2+...
              (zcoilbig{2}-rs(end,3)).^2);
%bigloop 5770
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment2(xcoilbig,ycoilbig,zcoilbig,numel(zcoilbig),[stind:numel(zcoilbig{numel(zcoilbig)})] ,h1,h2,rs,res,Fsym1...
    ,nrs);


Npts=200;
Npad=20;
stid=numel(rs(:,1))-160;
%rs=cat(1,rs,addnum2);
% 
% rs(stid-Npad:stid+Npts,1)= movmean(rs(stid-Npad:stid+Npts,1),Npad);
% rs(stid-Npad:stid+Npts,2)= movmean(rs(stid-Npad:stid+Npts,2),Npad);
% rs(stid-Npad:stid+Npts,3)= movmean(rs(stid-Npad:stid+Npts,3),Npad);
% rs(stid-Npad:stid+Npts,4)= movmean(rs(stid-Npad:stid+Npts,4),Npad);
% rs(stid-Npad:stid+Npts,5)= movmean(rs(stid-Npad:stid+Npts,5),Npad);
% 
% 
% Npts=200;
% Npad=40;
% stid=Npad+1;
% 
% rs(stid-Npad:stid+Npts,1)= movmean(rs(stid-Npad:stid+Npts,1),Npad);
% rs(stid-Npad:stid+Npts,2)= movmean(rs(stid-Npad:stid+Npts,2),Npad);
% rs(stid-Npad:stid+Npts,3)= movmean(rs(stid-Npad:stid+Npts,3),Npad);
% rs(stid-Npad:stid+Npts,4)= movmean(rs(stid-Npad:stid+Npts,4),Npad);
% rs(stid-Npad:stid+Npts,5)= movmean(rs(stid-Npad:stid+Npts,5),Npad);
close all
plot3(rs(:,1),rs(:,2),rs(:,3))
%pause
%%
close all
rsp=rs;

wireth1=.0021%+.0005;
wireth2=.0021%+.0005;
wirehe=.004;
close all
for i=1:numel(wid1)
subplot(3,4,i),plot3(rs(wid1(i):wid2(i),1),rs(wid1(i):wid2(i),2),rs(wid1(i):wid2(i),3))
axis equal
view(2)
end

%wireth=manualmodification(wireth1,wireth2,wid1,wid2,rs);
wireth=wireth1*ones(size(rs(:,1)));
rs=rsp;
np=numel(rs(:,1));
rs(end+1:2*end,1)=rs(1:np,1);
rs(np+1:2*np,3)=rs(1:np,3);
rs(np+1:2*np,4)=rs(1:np,4);
rs(np+1:2*np,2)=-rs(1:np,2);
rs(np+1:2*np,5)=-rs(1:np,5);
wireth(end+1:2*end)=wireth(1:np);

plot3(rs(:,1),rs(:,2),wireth);
axis equal


    len=sqrt((rs(2:end,1)-rs(1:end-1,1)).^2+...
             (rs(2:end,2)-rs(1:end-1,2)).^2+...
             (rs(2:end,3)-rs(1:end-1,3)).^2);
          [B,I]=sort(len);
          I=I(B==0);
          rs(I,:)=[];
          wireth(I)=[];
              len=sqrt((rs(2:end,1)-rs(1:end-1,1)).^2+...
             (rs(2:end,2)-rs(1:end-1,2)).^2+...
             (rs(2:end,3)-rs(1:end-1,3)).^2);
         len=cumsum([0;len(:)]);
    N=ceil(len(end)/.001);
    locs=0:len(end)/N:len(end);
    
     xcoil2=interp1(len,rs(1:end,1),locs,'line');
     ycoil2=interp1(len,rs(1:end,2),locs,'line');
     zcoil2=interp1(len,rs(1:end,3),locs,'line');
     wireth=interp1(len,wireth,locs,'line');
     
     xh=interp1(len,rs(1:end,4),locs,'line');
     yh=interp1(len,rs(1:end,5),locs,'line');
clear rs
 rs(:,1)=xcoil2;
 rs(:,2)=ycoil2;
 rs(:,3)=zcoil2;
save(strcat(outfile,'_sendtokop.mat'),'rs');
wireth=(wireth(2:end)+wireth(1:end-1)).'/2;
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
save(strcat(outfile,'_sendtokop.mat'),'rs');

[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);

save coilmesh.mat
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
% nl=numel(p(:,1));
% p2=p;
% p2(:,2)=-p2(:,2);
% te2p2=te2p+nl;
% p=cat(1,p,p2);
% te2p=cat(1,te2p,te2p2);
% reg=cat(1,reg(:),reg(:));
% portid=[1,1,2,2,1,1,2,2];
% port(end+1:end+4,:)=port(1:4,:)+nl;
% save coilmesh.mat te2p p reg port portid;
t2p=surftri(p,te2p);
np2=numel(p)/3;
p(end+1:end+8,1)=.01*[0,1,0,1,0,1,0,1];
p(np2+1:np2+8,2)=.01*[0,0,1,1,0,0,1,1];
p(np2+1:np2+8,3)=.01*[0,0,0,0,1,1,1,1]-.04;
t2p(end+1,:)=[np2+1;np2+2;np2+3];
t2p(end+1,:)=[np2+2;np2+4;np2+3];
t2p(end+1,:)=[np2+7;np2+6;np2+5];
t2p(end+1,:)=[np2+6;np2+7;np2+8];

t2p(end+1,:)=[np2+1;np2+2;np2+5];
t2p(end+1,:)=[np2+5;np2+6;np2+2];
t2p(end+1,:)=[np2+3;np2+4;np2+7];
t2p(end+1,:)=[np2+7;np2+8;np2+4];

t2p(end+1,:)=[np2+1;np2+3;np2+5];
t2p(end+1,:)=[np2+5;np2+7;np2+3];
t2p(end+1,:)=[np2+2;np2+4;np2+6];
t2p(end+1,:)=[np2+6;np2+8;np2+4];
load fourcorners.mat

t2p=cat(1,t2p,tri+numel(p(:,1)));
p=cat(1,p,pS);
trisurf(t2p,p(:,1),p(:,2),p(:,3));
tr=triangulation(t2p,p);
stlwrite(tr,strcat(outfile,'_coilwindingsrena.stl'))
addpath('C:\Users\ljg24\Desktop\optimisationcodes\anycoils\lfMFIE\wirecode')
addpath('C:\Users\ljg24\Desktop\optimisationcodes\anycoils\')
trisurf(tr);
%%

trisurf(t2p(end-11:end,:),p(:,1),p(:,2),p(:,3))
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
 

sigma=5.814*10^7;
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
jv=Nturns*Eout;

vol=mesh.tevol;
L=computeinductance(rv,jv)
p=mesh.p;
close all
 trisurf(tri,mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),'edgealpha',0,'facecolor',[127;127;127]/255)
 hold on
  trisurf(mesh.f2p(mesh.portf,:),mesh.p(:,1),mesh.p(:,2),mesh.p(:,3),double(mesh.portid),'edgealpha',0)
  axis equal
  caxis([1 2])
  light
  %%
[Jdc,rv,Ldc,R]=wirevieLaccurate_cl(meshfile,2*pi*3000,sigma,Nturns);
Jdc(:,1)=Jdc(:,1).*mesh.tevol(:);
Jdc(:,2)=Jdc(:,2).*mesh.tevol(:);
Jdc(:,3)=Jdc(:,3).*mesh.tevol(:);
   norm(Jdc(:)+jv(:))/norm(Jdc(:))
   
   
save(strcat(outfile,'_fdtmscoilsp.mat'),'rv','jv','tri','p','Jdc','R','Ldc');

end
 

function [xcoilb,ycoilb,zcoilb]=removenegside(xcoilbig,ycoilbig,zcoilbig)
keep=[];
for i=1:numel(ycoilbig)
if sum(ycoilbig{i})>0
    keep(end+1)=i;
end
end
for i=1:numel(keep)
xcoilb{i}=xcoilbig{keep(i)};
ycoilb{i}=ycoilbig{keep(i)};
zcoilb{i}=zcoilbig{keep(i)};
end
end



function [xcoil,ycoil,zcoil]=reorientcoil(xcoil,ycoil,zcoil);

for i=1:numel(xcoil)
            rct2(1)=mean(xcoil{i});
            rct2(2)=mean(ycoil{i});
            rct2(3)=mean(zcoil{i});
anglesx=atan2((ycoil{i}(:)-rct2(2)),(xcoil{i}(:)-rct2(1)));

anglesxtemp=-3;
    
[B,I]=sort(abs(anglesxtemp(1)-anglesx));
            if I(1)==1
    xcoil{i}=xcoil{i};
    ycoil{i}=ycoil{i};
    zcoil{i}=zcoil{i};
            else
    xcoil{i}=cat(2,xcoil{i}(I(1):end),xcoil{i}(1:I(1)-1));
    ycoil{i}=cat(2,ycoil{i}(I(1):end),ycoil{i}(1:I(1)-1));
    zcoil{i}=cat(2,zcoil{i}(I(1):end),zcoil{i}(1:I(1)-1));
            end

end
end



function [rs,nrs,stw,enw]=coilsegment(Vx,Vy,Vz,c1,c1ran,layer,layerp,rsp,res,Fsym,Vxt,Vyt,Vzt,c2,nrsp);
load firstordercorr_middle.mat nx ny nz
cent(1)=mean(Vxt{c2(1)});
cent(2)=mean(Vyt{c2(1)});
cent(3)=mean(Vzt{c2(1)});
cen(1)=mean(Vx{c1(1)});
cen(2)=mean(Vy{c1(1)});
cen(3)=mean(Vz{c1(1)});
cen=cen-cent;
len=cumsum(sqrt((Vx{c1(1)}(2:end)-Vx{c1(1)}(1:end-1)).^2+...
                (Vy{c1(1)}(2:end)-Vy{c1(1)}(1:end-1)).^2+...
                (Vz{c1(1)}(2:end)-Vz{c1(1)}(1:end-1)).^2));
len=[0;len(:)];
rang=min(len):(max(len)-min(len))/numel(Vxt{c1(1)}(:)):max(len);
[~,uind,~]=unique(len);

Vx{c1(1)}=interp1(len(uind),Vx{c1(1)}(uind),rang,'linear');
Vy{c1(1)}=interp1(len(uind),Vy{c1(1)}(uind),rang,'linear');
Vz{c1(1)}=interp1(len(uind),Vz{c1(1)}(uind),rang,'linear');
[~,I]=sort(sqrt((Vxt{c2(1)}(c1ran(1))-Vx{c1(1)}(:)+cen(1)).^2+...
                (Vyt{c2(1)}(c1ran(1))-Vy{c1(1)}(:)+cen(2)).^2+...
                (Vzt{c2(1)}(c1ran(1))-Vz{c1(1)}(:)+cen(3)).^2),1);
c1ran(1)
I(1)
            beg=I(1);
 
[~,I]=sort(sqrt((Vxt{c2(1)}(c1ran(end))-Vx{c1(1)}(:)+cen(1)).^2+...
                (Vyt{c2(1)}(c1ran(end))-Vy{c1(1)}(:)+cen(2)).^2+...
                (Vzt{c2(1)}(c1ran(end))-Vz{c1(1)}(:)+cen(3)).^2),1);
            
ending=I(1);
if ending<beg
c1ran=[beg:numel(Vx{c1(1)}) 1:ending];
else
c1ran=beg:ending;
end

load layerinfo.mat layerp;
rs(:,1)=Vx{c1(1)}(c1ran(:))+(layer)*nx(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,2)=Vy{c1(1)}(c1ran(:))+(layer)*ny(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,3)=Vz{c1(1)}(c1ran(:))+(layer)*nz(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,4)=Vx{c1(1)}(c1ran(:));
rs(:,5)=Vy{c1(1)}(c1ran(:));

nrs=floor(numel(rs)/8);

stw=1;
enw=numel(rs(:,1));
if numel(rsp)~=0

len=sqrt(sum((rsp(end,:)-rs(1,:)).^2,2));
N=ceil(len/res);
co=1/N:1/N:1-1/N;

rstemp(:,1)=(rsp(end,1)-layerp*nx(rsp(end,1),rsp(end,2)))*(1-co(:))+Vx{c1(1)}(c1ran(1))*(co(:));
rstemp(:,2)=(rsp(end,2)-layerp*ny(rsp(end,1),rsp(end,2)))*(1-co(:))+Vy{c1(1)}(c1ran(1))*(co(:));
rstemp(:,3)=Fsym(rstemp(:,1),rstemp(:,2));
rst(:,1)=rstemp(:,1)+(layerp*(1-co(:))+layer*co(:)).*nx(rstemp(:,1),rstemp(:,2));
rst(:,2)=rstemp(:,2)+(layerp*(1-co(:))+layer*co(:)).*ny(rstemp(:,1),rstemp(:,2));
rst(:,3)=rstemp(:,3)+(layerp*(1-co(:))+layer*co(:)).*nz(rstemp(:,1),rstemp(:,2));
rst(:,4)=rstemp(:,1);
rst(:,5)=rstemp(:,2);
Npts=numel(rst(:,1));
stid=numel(rsp(:,1));

 Npad=150;
if Npad>nrs || Npad>nrsp
     Npad=min(nrs,nrsp);
end
ending=stid+Npts+Npad;
rs=cat(1,rsp,rst,rs);
stw=numel(rsp(:,1))+numel(rst(:,1));
enw=numel(rs(:,1));
if stid-Npad>0
     rs(stid-Npad:ending,1)= movmean(rs(stid-Npad:ending,1),Npad);
 rs(stid-Npad:ending,2)= movmean(rs(stid-Npad:ending,2),Npad);
 rs(stid-Npad:ending,3)= movmean(rs(stid-Npad:ending,3),Npad);
 rs(stid-Npad:ending,4)= movmean(rs(stid-Npad:ending,4),Npad);
rs(stid-Npad:ending,5)= movmean(rs(stid-Npad:ending,5),Npad);
end
end
layerp=layer;
save layerinfo.mat layerp;
end

function [rs,nrs,stw,enw]=coilsegment2(Vx,Vy,Vz,c1,c1ran,layer,layerp,rsp,res,Fsym,nrsp);
load firstordercorr_middle.mat nx ny nz
load layerinfo.mat layerp;
rs(:,1)=Vx{c1(1)}(c1ran(:))+(layer)*nx(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,2)=Vy{c1(1)}(c1ran(:))+(layer)*ny(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,3)=Vz{c1(1)}(c1ran(:))+(layer)*nz(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,4)=Vx{c1(1)}(c1ran(:));
rs(:,5)=Vy{c1(1)}(c1ran(:));

nrs=floor(numel(rs)/8);

stw=1;
enw=numel(rs(:,1));
if numel(rsp)~=0

len=sqrt(sum((rsp(end,:)-rs(1,:)).^2,2));
N=ceil(len/res);
co=1/N:1/N:1-1/N;

rstemp(:,1)=(rsp(end,1)-layerp*nx(rsp(end,1),rsp(end,2)))*(1-co(:))+Vx{c1(1)}(c1ran(1))*(co(:));
rstemp(:,2)=(rsp(end,2)-layerp*ny(rsp(end,1),rsp(end,2)))*(1-co(:))+Vy{c1(1)}(c1ran(1))*(co(:));
rstemp(:,3)=Fsym(rstemp(:,1),rstemp(:,2));
rst(:,1)=rstemp(:,1)+(layerp*(1-co(:))+layer*co(:)).*nx(rstemp(:,1),rstemp(:,2));
rst(:,2)=rstemp(:,2)+(layerp*(1-co(:))+layer*co(:)).*ny(rstemp(:,1),rstemp(:,2));
rst(:,3)=rstemp(:,3)+(layerp*(1-co(:))+layer*co(:)).*nz(rstemp(:,1),rstemp(:,2));
rst(:,4)=rstemp(:,1);
rst(:,5)=rstemp(:,2);
Npts=numel(rst(:,1));
stid=numel(rsp(:,1));

 Npad=150;
if Npad>nrs || Npad>nrsp
     Npad=min(nrs,nrsp);
end
ending=stid+Npts+Npad;
rs=cat(1,rsp,rst,rs);
stw=numel(rsp(:,1))+numel(rst(:,1));
enw=numel(rs(:,1));
if stid-Npad>0
     rs(stid-Npad:ending,1)= movmean(rs(stid-Npad:ending,1),Npad);
 rs(stid-Npad:ending,2)= movmean(rs(stid-Npad:ending,2),Npad);
 rs(stid-Npad:ending,3)= movmean(rs(stid-Npad:ending,3),Npad);
 rs(stid-Npad:ending,4)= movmean(rs(stid-Npad:ending,4),Npad);
rs(stid-Npad:ending,5)= movmean(rs(stid-Npad:ending,5),Npad);
end
end
layerp=layer;
save layerinfo.mat layerp;
end