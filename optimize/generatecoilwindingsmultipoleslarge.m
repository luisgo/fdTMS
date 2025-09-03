function [xcoilbig,ycoilbig,zcoilbig,zzcoilbig,...
          xcoilfig8,ycoilfig8,zcoilfig8,zzcoilfig8,...
          xcoilside,ycoilside,zcoilside,zzcoilside]=...
          generatecoilwindingsmultipolesalllarge(harmfile,designfile,nlayers,desid,nloops,Ncc)
cutlev=[.07 -.07];
%generates a standard coil that is ultimately not used
[~,~,~,~,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,nlayers,desid,10,Ncc);
pause
%compute mesh normals
ntss=numel(t2psym)/12;
t2psym2=t2psym;
t2psym2(ntss+1:2*ntss,:)=t2psym(ntss+1:2*ntss,[1,3,2]);

t2psym2(2*ntss+1:3*ntss,:)=t2psym(2*ntss+1:3*ntss,[1,3,2]);
TR=triangulation(t2psym2,psym(:,1),psym(:,2),psym(:,3));
norp=vertexNormal(TR);


%plot for debugging
trisurf(TR,Att1,'edgealpha',0,'facecolor','interp')
%normalize stream function
normi=1/max(Att1(:));
Att1=Att1*normi;

%cut coil
[t2p,tind,Attmod,pmod,bdrmod]=TEMPgnerateholesnew(psym,t2psym,Att1,Fsym1,cutlev);

 close all

%get coil subregions
close all
nrr=1:numel(t2p);
for i=nrr
rr(i)=numel(t2p{i});
end
[val,~,rr2]=unique(rr);
reg{1}=nrr(rr2==3);
reg{2}=nrr(rr2==2);
reg{3}=nrr(rr2==1);

%combine common coil regions
for j=1:3
    t2p2{j}=t2p{reg{j}(1)};
    tind2{j}=tind{reg{j}(1)};
    pmod2{j}=pmod{reg{j}(1)};
    Attmod2{j}=Attmod{reg{j}(1)};
    bdrmod2{j}=bdrmod{reg{j}(1)};
    for i=2:numel(reg{j})
        bdrmod2{j}=cat(1, bdrmod2{j}(:),bdrmod{reg{j}(i)}(:)+max(t2p2{j}(:)));
        t2p2{j}=cat(1,t2p2{j},t2p{reg{j}(i)}+max(t2p2{j}(:)));
        tind2{j}=cat(1,tind2{j},tind{reg{j}(i)});
        pmod2{j}=cat(1, pmod2{j},pmod{reg{j}(i)});
        Attmod2{j}=cat(1, Attmod2{j}(:),Attmod{reg{j}(i)}(:));
    end
end
t2p=t2p2; 

%makesure figure8 region is 2
if max(abs(Attmod2{3}))>.7 
t2ptemp=t2p{2};
t2p{2}=t2p{3};
t2p{3}=t2ptemp;

t2ptemp=tind2{2};
tind2{2}=tind2{3};
tind2{3}=t2ptemp;

t2ptemp=pmod2{2};
pmod2{2}=pmod2{3};
pmod2{3}=t2ptemp;


t2ptemp=Attmod2{2};
Attmod2{2}=Attmod2{3};
Attmod2{3}=t2ptemp;


t2ptemp=bdrmod2{2};
bdrmod2{2}=bdrmod2{3};
bdrmod2{3}=t2ptemp;

end
for i=1:numel(t2p)
hold on
trisurf(t2p{i},pmod2{i}(:,1),pmod2{i}(:,2),pmod2{i}(:,3),Attmod2{i}(:),'edgealpha',1,'facealpha',1,'facecolor','interp')
end

pause
t2pmod2=t2p;


%%%%%%start of actual coil design
load(harmfile)
load(designfile);
%%%%make normal

nx = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,1));
ny = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,2));
nz = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,3));
%%%where to evaluate field
load barisamp4.mat
addpath('C:\Users\ljg24\Desktop\FEM_modes\matlab');
p=.073*p;
        js(:,3)=Az(:,1:Ncc)*Jcalc2{desid};
        js(:,2)=Ay(:,1:Ncc)*Jcalc2{desid};
        js(:,1)=Ax(:,1:Ncc)*Jcalc2{desid};
%close all
%quiver3(rs(:,1),rs(:,2),rs(:,3),js(:,1),js(:,2),js(:,3))
[Eoutreg]=computeEprimary(loc',js',numel(js)/3,p',numel(p)/3);
max(abs(Eoutreg(:)))
levs=-1+1/nloops:2/nloops:1-1/nloops
normi
[rs,js,Eout,p,~,~,...
    ~,~,~,~]=coilfields(psym,t2psym,Att1,levs);

mean(abs(Eoutreg(:)))/mean(abs(Eout(:)))

normi=-mean(abs(Eoutreg(:)))/mean(abs(Eout(:)))/2;
max(Eout(:))
norm(Eout(:)-Eoutreg(:))/norm(Eoutreg(:))
addpath('C:\workbackup\ljg24\loopstar\Arraycode\Simtools\InductanceEnergy')

%%%%big loop optimization
region=1;
Nv=[3;4;1];
Aine=[];
X0=[];
LB=[];
UB=[];
i=2;
for i=1:3
nhhmod{i}(:,1)=nx(pmod2{i}(:,1),pmod2{i}(:,2));
nhhmod{i}(:,2)=ny(pmod2{i}(:,1),pmod2{i}(:,2));
nhhmod{i}(:,3)=nz(pmod2{i}(:,1),pmod2{i}(:,2));
maii(i)=max(abs(Attmod2{i}))-.005;
mii(i)=getmin2(pmod2{i},t2pmod2{i},Attmod2{i},i)
X0p=linspace(mii(i),maii(i),Nv(i));
A=eye(Nv(i)-1);
cdist{i}=getmin(pmod2{i},t2pmod2{i},Attmod2{i},Fsym1,mii(i),maii(i));
for j=1:Nv(i)-1;
A(j,j+1)=-1;
end
size(A)
Aine(end+1:end+Nv(i)-1,end+1:end+Nv(i))=A;
X0=cat(1,X0(:),X0p(:));
LB=cat(1,LB(:),mii(i)*ones(Nv(i),1));
UB=cat(1,UB(:),maii(i)*ones(Nv(i),1));
end
Bine=zeros([sum(Nv)-numel(Nv),1]);
 Aine(:,end+1)=0;
 LB(end+1)=normi*(1+.5);
 UB(end+1)=normi*(1-.5);
 X0(end+1)=normi;  
%Jval=sqrt(2*Menergy(desid)/(10^-5));
options = optimoptions('fmincon');
options.MaxIterations=100;
options.ObjectiveLimit=0.03;
X = fmincon(@(levs)costfnalllarge(levs(1:end-1),Eoutreg,pmod2,t2pmod2,Attmod2,levs(end),Nv,nhhmod),X0,Aine,Bine,[],[],...
    LB,UB,@(levs) mycon(levs(1:end-1),cdist,Nv),options);
normi=X(end)
c=costfnalllarge(X(1:end-1),Eoutreg,pmod2,t2pmod2,Attmod2,normi,Nv,nhhmod)
Z=X(Nv(1)+1:Nv(1)+Nv(2));
Y=X(Nv(1)+Nv(2)+1:sum(Nv));
X=X(1:Nv(1));

levs=-1+1/nloops:2/nloops:1-1/nloops
for i=1:3
if i==1
[rs1,js1,Eout,p,t2p,clen,...
    xcoilbig,ycoilbig,zcoilbig,zzcoilbig]=coilfields(pmod2{i},t2pmod2{i},Attmod2{i},cat(1,X(:),-X(:)));
elseif  i==2
[rs2,js2,Eout,p,t2p,clen,...
    xcoilfig8,ycoilfig8,zcoilfig8,zzcoilfig8]=coilfields(pmod2{i},t2pmod2{i},Attmod2{i},cat(1,Z(:),-Z(:)));
elseif  i==3
[rs3,js3,~,p,t2p,clen,...
    xcoilside,ycoilside,zcoilside,zzcoilside]=coilfields(pmod2{i},t2pmod2{i},Attmod2{i},cat(1,Y(:),-Y(:)));
end
end
end

function cc=getmin(p,t2p,Att,Fsym1,Imin,Imax)

dfun=Imin:(Imax-Imin)/50:Imax;
[dfun1 dfun2]=ndgrid(dfun,dfun);
yfun=1000*ones([numel(dfun) numel(dfun)]);
for i=1:numel(dfun)
    for j=1:numel(dfun)
[xcoil ycoil zcoil]=getloop([dfun(i) dfun(i)],p,t2p,Att,Fsym1,0);
[xcoil2 ycoil2 zcoil2]=getloop([dfun(j) dfun(j)],p,t2p,Att,Fsym1,1);

for k=1:numel(xcoil)
ddd=min((xcoil(k)-xcoil2(:)).^2+...
        (ycoil(k)-ycoil2(:)).^2+...
        (zcoil(k)-zcoil2(:)).^2);
    if ddd<yfun(i,j)
        yfun(i,j)=ddd;
    end
end

    end
end

cc=scatteredInterpolant(dfun1(:),dfun2(:),sqrt(yfun(:)),'natural'); 
end

function  [c,ceq]=mycon(levs,cdist,Nv)
c=zeros([sum(Nv)-numel(Nv) 1]);
ct=1;
for j=1:numel(Nv);
    stn=sum(Nv(1:j-1));
for i=1:Nv(j)-1
c(ct)=-cdist{j}(levs(stn+i),levs(stn+i+1))+0.003;
ct=ct+1;
end
end
ceq=[];
end

function [xcoil,ycoil,zcoil]=getloop(lev,p,t2p,Att,Fsym1,inte)
c = tricontour(p(:,1:2),t2p,Att,[lev lev]);

ct=1;
en=0;
xcoil=[];ycoil=[];zcoil=[];

while en<numel(c(:))/2
st=en+1;
en=c(2,st)+st;
height=c(1,st);
st=st+1;

if inte==1
c3=Fsym1(c(1,st:en),c(2,st:en));
    len=sqrt((c(1,st+1:en)-c(1,st:en-1)).^2+...
             (c(2,st+1:en)-c(2,st:en-1)).^2+...
             (c3(1,2:end)-c3(1,1:end-1)).^2);         
     len=cumsum([0;len(:)])';
     N=ceil(len(end)/.0001);
     locs=0:len(end)/N:len(end);
     xcoil=cat(1,xcoil,interp1(len,c(1,st:en),locs,'linear')');
     ycoil=cat(1,ycoil,interp1(len,c(2,st:en),locs,'linear')');
     zcoil=cat(1,zcoil,interp1(len,c3(1,1:end),locs,'linear')');
     else
xcoil=cat(1,xcoil,c(1,st:en)');
ycoil=cat(1,ycoil,c(2,st:en)');
zcoil=cat(1,zcoil,Fsym1(c(1,st:en),c(2,st:en))');
end
ct=ct+1;
end
end
    


function cc=getmin2(p,t2p,Att,ival)
pcen=(p(t2p(:,1),:)+p(t2p(:,2),:)+p(t2p(:,3),:))/3;
if ival<3
lo=(pcen(:,2)>0);
else
lo=(pcen(:,1)>0).*(pcen(:,2)>0);    
end
t2p=t2p(lo==1,:);
bat=min(abs(Att)):.001:min(abs(Att))+.2;
cc=10;
for i=bat
ct=getloop2(i,p,t2p,abs(Att));
if ct==1 
if i<cc
cc=i;
break;
end
end
end

bat=cc-.001:.00001:cc;

for i=bat
ct=getloop2(i,p,t2p,abs(Att));
if ct==1
if i<cc
cc=i;
break;
end
end
end

end


    function ct=getloop2(lev,p,t2p,Att)
cou = tricontour(p(:,1:2),t2p,Att,[lev lev]);
ct=0;
en=0;
c=cou;
if numel(cou)==0
ct=2;
return;
end
while en<numel(c(1,:))
st=en+1;
en=c(2,st)+st;
st=st+1;
ct=ct+1;
end
if ct==1
if norm(c(:,st)-c(:,en))>.001
    ct=2;
end
if min(abs(c(2,st:en)))<=.0015
    ct=2;
end
end
    end

