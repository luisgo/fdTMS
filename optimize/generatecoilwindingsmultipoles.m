function [xcoilbig,ycoilbig,zcoilbig,zzcoilbig,...
          xcoilfig8,ycoilfig8,zcoilfig8,zzcoilfig8,...
          xcoilside,ycoilside,zcoilside,zzcoilside]=...
          generatecoilwindingsmultipoles(harmfile,designfile,nlayers,desid,nloops,Ncc)

%generates a standard coil that is ultimately not used
[~,~,~,~,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,nlayers,desid,10,Ncc);
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
Att1=Att1/max(Att1);
normi=1/max(Att1(:));
%cut coil
[t2p,tind]=TEMPgnerateholes(psym,t2psym,Att1,Fsym1,[.07 -.07]);
[t2p,tind,Attmod,pmod]=TEMPgnerateholesnew(psym,t2psym,Att1,Fsym1,[.07 -.07]);

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
    for i=2:numel(reg{j})
        t2p2{j}=cat(1,t2p2{j},t2p{reg{j}(i)});
        tind2{j}=cat(1,tind2{j},tind{reg{j}(i)});
        pmod2{j}(unique(t2p{reg{j}(i)}),:)=pmod{reg{j}(i)}(unique(t2p{reg{j}(i)}),:);
    end
end
t2p=t2p2; 
pss=unique(t2p{3});
%makesure figure8 region is 2
if max(abs(Att1(pss)))>.7 
t2ptemp=t2p{2};
t2p{2}=t2p{3};
t2p{3}=t2ptemp;
t2ptemp=tind2{2};
tind2{2}=tind2{3};
tind2{3}=t2ptemp;
end
for i=1:numel(t2p)
 %subplot(3,3,i),
hold on
col=i*ones([numel(t2p{i})/3 1]); 
trisurf(t2p{i},pmod2{i}(:,1),pmod2{i}(:,2),pmod2{i}(:,3),col,'edgealpha',0,'facealpha',.4)
text(psym(t2p{i}(1,1),1),psym(t2p{i}(1,1),2),psym(t2p{i}(1,1),3),num2str(i));
end
pause
t2preg=t2p;
%%%%%%start of actual coil design
load(harmfile)
load(designfile);
j=1;
%%%%make normal

nx = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,1));
ny = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,2));
nz = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,3));
%%%where to evaluate field
load barisamp4.mat
addpath('C:\Users\ljg24\Desktop\FEM_modes\matlab');
p=.073*p;

%%%%% get fields for region 1
levs=-1+1/nloops:2/nloops:1-1/nloops
[rs,js]=coilfields(psym,t2preg{1},Att1,levs);
nhh=zeros(size(rs));
nhh(:,1)=nx(rs(:,1),rs(:,2));
nhh(:,2)=ny(rs(:,1),rs(:,2));
nhh(:,3)=nz(rs(:,1),rs(:,2));
rs=cat(1,rs-.002*nhh,rs+.002*nhh);
js=cat(1,js,js);
%js=2*js;
[Eoutreg1]=computeEprimary(rs',js',numel(rs)/3,p',numel(p)/3);

%%%%% get fields for region 2
levs=-1+1/nloops:2/nloops:1-1/nloops
[rs,js]=coilfields(psym,t2preg{2},Att1,levs);
nhh=zeros(size(rs));
nhh(:,1)=nx(rs(:,1),rs(:,2));
nhh(:,2)=ny(rs(:,1),rs(:,2));
nhh(:,3)=nz(rs(:,1),rs(:,2));
rs=cat(1,rs-.002*nhh,rs+.002*nhh);
js=cat(1,js,js);
%js=2*js;
[Eoutreg2]=computeEprimary(rs',js',numel(rs)/3,p',numel(p)/3);

%%%%% get fields for region 3
[rs,js]=coilfields(psym,t2preg{3},Att1,levs);
nhh=zeros(size(rs));
nhh(:,1)=nx(rs(:,1),rs(:,2));
nhh(:,2)=ny(rs(:,1),rs(:,2));
nhh(:,3)=nz(rs(:,1),rs(:,2));
rs=cat(1,rs-.002*nhh,rs+.002*nhh);
js=cat(1,js,js);
%js=2*js;
[Eoutreg3]=computeEprimary(rs',js',numel(rs)/3,p',numel(p)/3);

addpath('C:\workbackup\ljg24\loopstar\Arraycode\Simtools\InductanceEnergy')

%%%%big loop optimization
region=1;
Nv=3;
A=eye(Nv-1);
for i=1:Nv-1;
A(i,i+1)=-1;
end
B=zeros([Nv-1,1]);
mii=min(abs(Att1(t2preg{region}(:))))+.006;
maii=max(abs(Att1(t2preg{region}(:))))-.01;
% [t2preg{region},tind2{region},nodes]=padtriag(psym,t2psym,t2preg{region});
% [t2preg{region},tind2{region},nodes]=padtriag(psym,t2psym,t2preg{region});
% tind2{region}=setdiff(tind2{region},tind2{2});
% tind2{region}=setdiff(tind2{region},tind2{3});
% t2preg{region}=t2psym(tind2{region},:);
X0=linspace(mii,maii,Nv);
psym2=psym-.002*norp;
X = fmincon(@(levs)costfn(levs,Eoutreg1,psym2,t2preg{region},Att1),X0,A,B,[],[],...
    mii*ones(Nv,1),maii*ones(Nv,1));
[rs,js,Eout1,p,t2p,clen]=coilfields(psym2,t2preg{region},Att1,cat(1,X(:),-X(:)));
 Ldc=computeinductance(rs,js)
 Rdc=clen/(4*2*58.14)

 
 
 
%%%%fig8 loop optimization

nx = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,1));
ny = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,2));
nz = scatteredInterpolant(psym(:,1),psym(:,2),norp(:,3));
region=2;
Nv=5;
A=eye(Nv-1);
for i=1:Nv-1;
A(i,i+1)=-1;
end
B=zeros([Nv-1,1]);
mii=min(abs(Att1(t2preg{region}(:))))+.004;
maii=max(abs(Att1(t2preg{region}(:))))-.004;
 [t2preg{region},tind2{region},nodes]=padtriag(psym,t2psym,t2preg{region});
 [t2preg{region},tind2{region},nodes]=padtriag(psym,t2psym,t2preg{region});
 tind2{region}=setdiff(tind2{region},tind2{1});
 tind2{region}=setdiff(tind2{region},tind2{3});
 t2preg{region}=t2psym(tind2{region},:);

 X0=linspace(mii,maii,Nv);
Z = fmincon(@(levs)costfn2(levs,Eoutreg2,psym,t2preg{region},Att1,nx,ny,nz),X0,A,B,[],[],...
    mii*ones(Nv,1),maii*ones(Nv,1));
[rs,js,Eout2,p,t2p,clen]=coilfieldsp(psym,t2preg{region},Att1,cat(1,Z(:),-Z(:)),nx,ny,nz);
 
%%%%side loop optimization
region=3;
Nv=2;
A=eye(Nv-1);
for i=1:Nv-1;
A(i,i+1)=-1;
end
B=zeros([Nv-1,1]);
maii=max(abs(Att1(t2preg{region}(:))))-.004;
[t2preg{region},tind2{region},nodes]=padtriag(psym,t2psym,t2preg{region});
[t2preg{region},tind2{region},nodes]=padtriag(psym,t2psym,t2preg{region});
tind2{region}=setdiff(tind2{region},tind2{2});
t2preg{region}=t2psym(tind2{region},:);
mii=getmin(psym,t2preg{region},Att1)+.001;

X0=linspace(mii+.05,maii,Nv);
psym2=psym-.002*norp;
Y = fmincon(@(levs)costfn(levs,Eoutreg3,psym2,t2preg{region},Att1),X0,A,B,[],[],...
    mii*ones(Nv,1),maii*ones(Nv,1));
[rs,js,Eout3,p,t2p,clen]=coilfields(psym2,t2preg{region},Att1,cat(1,Y(:),-Y(:)));



 
 
err1=norm(Eoutreg1(:)-Eout1(:))/norm(Eoutreg1(:))
err2=norm(Eoutreg2(:)-Eout2(:))/norm(Eoutreg2(:))
err3=norm(Eoutreg3(:)-Eout3(:))/norm(Eoutreg3(:))


levs=-1+1/nloops:2/nloops:1-1/nloops
for i=1:3
if i==1
[rs1,js1,Eout,p,t2p,clen,...
    xcoilbig,ycoilbig,zcoilbig,zzcoilbig]=coilfields(psym,t2preg{i},Att1,cat(1,X(:),-X(:)));
elseif  i==2
[rs2,js2,Eout,p,t2p,clen,...
    xcoilfig8,ycoilfig8,zcoilfig8,zzcoilfig8]=coilfields(psym,t2preg{i},Att1,levs);
elseif  i==3
[rs3,js3,Eout,p,t2p,clen,...
    xcoilside,ycoilside,zcoilside,zzcoilside]=coilfields(psym,t2preg{i},Att1,cat(1,Y(:),-Y(:)));
end
end
end


function cc=getmin(p,t2p,Att)
pcen=(p(t2p(:,1),:)+p(t2p(:,2),:)+p(t2p(:,3),:))/3;
lo=(pcen(:,1)>0).*(pcen(:,2)>0);
t2p=t2p(lo==1,:);
bat=min(abs(Att)):.001:min(abs(Att))+.1;
cc=10;
for i=bat
ct=getloop(i,p,t2p,abs(Att));
if ct==1
if i<cc
cc=i;
break;
end
end
end

bat=cc-.001:.00001:cc;

for i=bat
ct=getloop(i,p,t2p,abs(Att));
if ct==1
if i<cc
cc=i;
break;
end
end
end

end


    function ct=getloop(lev,p,t2p,Att)
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
ct=ct+1;
end
    end

