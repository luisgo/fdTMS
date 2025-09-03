load('hhmeshingkit.mat')
filename='renamesh';

for i=1:20
plot3(rs(wid1(i):wid2(i),1),rs(wid1(i):wid2(i),2),rs(wid1(i):wid2(i),3))
hold on
text(rs(wid1(i),1),rs(wid1(i),2),rs(wid1(i),3),num2str(i))
end

%find mesh dimensions
wireth1=.00675+.008%+.0005;
wireth2=.00675+.008%+.0005;
nor1=.007%+.0005;
nor2=.015%+.0005;
wirethick=manualmodification3layers(wireth1,wireth2,wid1,wid2,rs);
wirenorfig8=manualmodification3layers(nor2,nor2,wid1,wid2,rs);
wirenorside=wirenorfig8;

    wirenorside(wid1(1):wid2(4))=nor1;
    wirenorside(wid1(6):wid2(7))=nor1;
    wirenorside(wid1(10):wid2(9))=nor1;
    wirenorside(wid1(20):wid2(20))=nor1;


wireth1=.003%+.0005;
wireth2=.00675%+.0005;
wireth=manualmodification3layers(wireth1,wireth2,wid1,wid2,rs);
%mirror the mesh
rsp=rs;
rs=rsp;
np=numel(rs(:,1));
rs(end+1:2*end,1)=rs(1:np,1);
rs(np+1:2*np,3)=rs(1:np,3);
rs(np+1:2*np,4)=rs(1:np,4);
rs(np+1:2*np,2)=-rs(1:np,2);
rs(np+1:2*np,5)=-rs(1:np,5);
wirethick(end+1:2*end)=wirethick(1:np);
wirenorfig8(end+1:2*end)=wirenorfig8(1:np);
wirenorside(end+1:2*end)=wirenorside(1:np);
wireth(end+1:2*end)=wireth(1:np);

%% resample rs


    len=sqrt((rs(2:end,1)-rs(1:end-1,1)).^2+...
             (rs(2:end,2)-rs(1:end-1,2)).^2+...
             (rs(2:end,3)-rs(1:end-1,3)).^2);
          [B,I]=sort(len);
          I=I(B==0);
          rs(I,:)=[];
          wireth(I)=[];
     wirethick(I)=[];
     wirenorfig8(I)=[];
     wirenorside(I)=[];
              len=sqrt((rs(2:end,1)-rs(1:end-1,1)).^2+...
             (rs(2:end,2)-rs(1:end-1,2)).^2+...
             (rs(2:end,3)-rs(1:end-1,3)).^2);
         len=cumsum([0;len(:)]);
    N=ceil(len(end)/.002);
    locs=0:len(end)/N:len(end);
    
     xcoil2=interp1(len,rs(1:end,1),locs,'line');
     ycoil2=interp1(len,rs(1:end,2),locs,'line');
     zcoil2=interp1(len,rs(1:end,3),locs,'line');
     wireth=interp1(len,wireth,locs,'line');
     wirethick=interp1(len,wirethick,locs,'line');
     wirenorfig8=interp1(len,wirenorfig8,locs,'line');
     wirenorside=interp1(len,wirenorside,locs,'line');
     
     xh=interp1(len,rs(1:end,4),locs,'line');
     yh=interp1(len,rs(1:end,5),locs,'line');

     %% generate meshes
     

     clear rs
 rs(:,1)=xcoil2;
 rs(:,2)=ycoil2;
 rs(:,3)=zcoil2;
 
wireth=(wireth(2:end)+wireth(1:end-1)).'/2;
wirethick=(wirethick(2:end)+wirethick(1:end-1)).'/2;
wirenorfig8=(wirenorfig8(2:end)+wirenorfig8(1:end-1)).'/2;
wirenorside=(wirenorside(2:end)+wirenorside(1:end-1)).'/2;

drs=rs(2:end,:)-rs(1:end-1,:);
rs1=(rs(2:end,:)+rs(1:end-1,:))/2;
rs=rs1;
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

%thick wire mesh
[rs1,rs2,rs3,rs4]=givethickness(wirethick,wirenorfig8,tan,nhat,rs);
[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);
close all

t2p=surftri(p,te2p);
[p,t2p]=addcube(p,t2p,-0.04);
trisurf(t2p,p(:,1),p(:,2),p(:,3));
tr=triangulation(t2p,p);
  stlwrite(tr,strcat(filename,'thickwire.stl'))

pause
%ribbed wire mesh

[rs1,rs2,rs3,rs4]=givethickness(wireth,wirenorfig8,tan,nhat,rs);
[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);
close all
t2p=surftri(p,te2p);

t2p=surftri(p,te2p);
load C:\Users\ljg24\Desktop\shamcoil\generatecoilmesh\fourcorners.mat

t2p=cat(1,t2p,tri+numel(p(:,1)));
p=cat(1,p,pS);

[p,t2p]=addcube(p,t2p,-0.04);
trisurf(t2p,p(:,1),p(:,2),p(:,3));
tr=triangulation(t2p,p);
  stlwrite(tr,strcat(filename,'ribbedwire.stl'))
pause
  %stepped wire mesh

[rs1,rs2,rs3,rs4]=givethickness(wirethick,wirenorside,tan,nhat,rs);
[te2p,p,reg,port,portid,hexmesh]=maketetramesh(rs1,rs2,rs3,rs4);

t2p=surftri(p,te2p);
[p,t2p]=addcube(p,t2p,-0.04);
trisurf(t2p,p(:,1),p(:,2),p(:,3));
tr=triangulation(t2p,p);
  stlwrite(tr,strcat(filename,'steppedwire.stl'))

axis equal

