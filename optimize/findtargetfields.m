function [Eout,Eout1,Eout2,Eout3,Eoutr]=...
          findtargetfields(harmfile,designfile,nlayers,desid,nloops,Ncc)

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
[t2p,tind]=TEMPgnerateholes(psym,t2psym,Att1,Fsym1,[.07 -.07])
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
    for i=2:numel(reg{j})
        t2p2{j}=cat(1,t2p2{j},t2p{reg{j}(i)});
        tind2{j}=cat(1,tind2{j},tind{reg{j}(i)});

    end
end
t2p=t2p2; 


%%%%%%actualfields
load(harmfile)
load(designfile);
j=1;
Js(:,1)=Ax(:,1:Ncc)*Jcalc2{desid};
Js(:,2)=Ay(:,1:Ncc)*Jcalc2{desid};
Js(:,3)=Az(:,1:Ncc)*Jcalc2{desid};

load barisamp4.mat
addpath('C:\Users\ljg24\Desktop\FEM_modes\matlab');
p=.07*p;
Js=2*pi*3000*Js;

[Eout]=computeEprimary(loc(:,:)',Js(:,:)',numel(loc(:,:))/3,p',numel(p)/3);


levs=-1+1/nloops:2/nloops:1-1/nloops
[rs,js,Eoutr]=coilfields(psym,t2psym,Att1,levs);
[rs,js,Eout1]=coilfields(psym,t2p2{1},Att1,levs);
[rs,js,Eout2]=coilfields(psym,t2p2{2},Att1,levs);
[rs,js,Eout3]=coilfields(psym,t2p2{3},Att1,levs);
Eoutr=Eoutr*normi;
Eout1=Eout1*normi;
Eout2=Eout2*normi;
Eout3=Eout3*normi;

%[Eout1]=computeEprimary(rs',js',numel(loc(:,:))/3,p',numel(p)/3);




[Eout1]=computeEprimary(loc(tind2{1},:)',Js(tind2{1},:)',numel(loc(tind2{1},:))/3,p',numel(p)/3);
[Eout2]=computeEprimary(loc(tind2{2},:)',Js(tind2{2},:)',numel(loc(tind2{2},:))/3,p',numel(p)/3);
[Eout3]=computeEprimary(loc(tind2{3},:)',Js(tind2{3},:)',numel(loc(tind2{3},:))/3,p',numel(p)/3);
