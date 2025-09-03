clear all
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

%%inner middle coil
h1=-.004067;
h2=-.004067;
c1=[1,2,3,5,7];
%side L
c2=[8];
%side R
c3=[9];
%bottom
c3=[4,6];
res=.0001


 wid1=[];
 wid2=[];
layerp=h1;
save layerinfo.mat layerp;
nrs=100000;
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,6,[5800:numel(xcoil{6}) 1:1473] ,h1,h1,[],res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
%feed 1


Npts=50;
wid2(1)=wid2(1)+50;
 convi=(1/Npts:1/Npts:1-1/Npts);
 Npts=Npts-1;
 ax=rs(1,1)*ones(size(convi));ay=rs(1,2)*(convi);
 az=rs(1,3)*ones(size(convi));
  load ../modeheight.mat
[TH,PHI,R] = cart2sph(p(:,1),-p(:,2),p(:,3));
rad=scatteredInterpolant(TH,PHI,R);
[TH,PHI] = cart2sph(ax,ay,az);
[ax,ay,az]=sph2cart(TH,PHI,rad(TH,PHI));%project
 addnum(1:Npts,1)=ax+h1*nx(ax,ay);
 addnum(1:Npts,2)=ay+h1*ny(ax,ay);
 %addnum(1:Npts,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum(1:Npts,3)=az+h1*nz(ax,ay);
 addnum(1:Npts,4)=ax;
 addnum(1:Npts,5)=ay;
 stid=5;

  addnum2(1:Npts,1)=ax+h2*nx(ax,ay);
 addnum2(1:Npts,2)=ay+h2*ny(ax,ay);
 %addnum(1:Npts,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum2(1:Npts,3)=az+h2*nz(ax,ay);
 addnum2(1:Npts,4)=ax;
 addnum2(1:Npts,5)=ay;
 addnum2=addnum2(end:-1:1,:);

 plot3(addnum(:,1),addnum(:,2),addnum(:,3))
 hold on
 
 plot3(addnum2(:,1),addnum2(:,2),addnum2(:,3),'--r')
 axis equal
 view(2)
 rs=cat(1,addnum(1:end-5,:),rs(5:end,:));
 % %round corners
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,4,[673:numel(xcoil{4}) 1:638] ,h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
[rs,nrs]=coilsegment(xcoil,ycoil,zcoil,4,[628:numel(xcoil{4}) 1:590] ,h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,6,[1564:4566],h1,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
%y=29 to 32
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,7,[455:703],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,9,[221:numel(xcoil{9}) 1:171] ,h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
 [rs,nrs]=coilsegment(xcoil,ycoil,zcoil,9,[215:numel(xcoil{9}) 1:120] ,h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
 [rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,7,[783:numel(xcoil{7}) 1:100],h1,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);

 
 [rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,8,[579:numel(xcoil{8}) 1:509],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
[rs,nrs]=coilsegment(xcoil,ycoil,zcoil,8,[535:numel(xcoil{8}) 1:448],h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
 [rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,7,[214:419],h1,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
% 


[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,5,[378:numel(xcoil{5}) 1:337],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,3,[261:numel(xcoil{3}) 1:220],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,2,[230:numel(xcoil{2}) 1:189],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);


%[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,1,[115:numel(xcoil{1}) 1:74],h1,h1,rs,res,Fsym1...
%    ,xcoil2,ycoil2,zcoil2,nrs);
 addnum=rs(end,:);
 convi=(1/30:1/30:1);

  ax=rs(end,1)*(1-convi);ay=rs(end,2)-.005*convi;
 addnum(1:30,1)=ax+h1*nx(ax,ay);
 addnum(1:30,2)=ay+h1*ny(ax,ay);
 addnum(1:30,3)=Fsym1(ax,ay)+h1*nz(ax,ay);
 addnum(1:30,4)=ax;
 addnum(1:30,5)=ay;
 
 rs=cat(1,rs,addnum);
 
 [rs,nrs]=coilsegment(xcoil,ycoil,zcoil,1,[115:numel(xcoil{1}) 1:95],h1,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);

 [rs,nrs]=coilsegment(xcoil,ycoil,zcoil,1,[115:numel(xcoil{1}) 1:74],h2,h1,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
 [rs,nrs]=coilsegment(xcoil,ycoil,zcoil,2,[230:numel(xcoil{2}) 1:189],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
 [rs,nrs]=coilsegment(xcoil,ycoil,zcoil,3,[261:numel(xcoil{3}) 1:220],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
 [rs,nrs]=coilsegment(xcoil,ycoil,zcoil,5,[378:numel(xcoil{5}) 1:337],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
 [rs,nrs]=coilsegment(xcoil,ycoil,zcoil,7,[460:numel(xcoil{7}) 1:429],h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
[rs,nrs,wid1(end+1),wid2(end+1)]=coilsegment(xcoil,ycoil,zcoil,6,[4616:5770] ,h1,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);
[rs,nrs]=coilsegment(xcoil,ycoil,zcoil,6,[5800:numel(xcoil{6}) 1:5750] ,h2,h2,rs,res,Fsym1...
    ,xcoil2,ycoil2,zcoil2,nrs);

Npts=200;
Npad=20;
stid=numel(rs(:,1))-160;
 rs=cat(1,rs,addnum2);

rs(stid-Npad:stid+Npts,1)= movmean(rs(stid-Npad:stid+Npts,1),Npad);
rs(stid-Npad:stid+Npts,2)= movmean(rs(stid-Npad:stid+Npts,2),Npad);
rs(stid-Npad:stid+Npts,3)= movmean(rs(stid-Npad:stid+Npts,3),Npad);
rs(stid-Npad:stid+Npts,4)= movmean(rs(stid-Npad:stid+Npts,4),Npad);
rs(stid-Npad:stid+Npts,5)= movmean(rs(stid-Npad:stid+Npts,5),Npad);


Npts=200;
Npad=40;
stid=Npad+1;

rs(stid-Npad:stid+Npts,1)= movmean(rs(stid-Npad:stid+Npts,1),Npad);
rs(stid-Npad:stid+Npts,2)= movmean(rs(stid-Npad:stid+Npts,2),Npad);
rs(stid-Npad:stid+Npts,3)= movmean(rs(stid-Npad:stid+Npts,3),Npad);
rs(stid-Npad:stid+Npts,4)= movmean(rs(stid-Npad:stid+Npts,4),Npad);
rs(stid-Npad:stid+Npts,5)= movmean(rs(stid-Npad:stid+Npts,5),Npad);

plot3(rs(:,1),rs(:,2),rs(:,3));
hold on
trisurf(t2p,p(:,1),-p(:,2),p(:,3),'edgealpha',0,'facealpha',.5)
rsp=rs;


wireth1=.00575+.008%+.0005;
wireth2=.00575+.008%+.0005;
wirehe=.02467;
close all
for i=1:numel(wid1)
subplot(3,4,i),plot3(rs(wid1(i):wid2(i),1),rs(wid1(i):wid2(i),2),rs(wid1(i):wid2(i),3))
axis equal
view(2)
end

wireth=manualmodification(wireth1,wireth2,wid1,wid2,rs);
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
    N=ceil(len(end)/.003);
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
save sendtokop.mat rs;
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
save sendtokop.mat rs;

[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);


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
p(np2+1:np2+8,3)=.01*[0,0,0,0,1,1,1,1];
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
stlwrite(tr,'..\coilwindingsrena2.stl')
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