clear all
harmfile='../newwireharmonics.mat';
dept=142;Wc=107;Sc=15.2;
Ncc=16;
designfile=strcat('optimT3d',num2str(dept),'ncoils16lmax30.mat');
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    ct=1;
i=3;
[xcoilc{ct},ycoilc{ct},zcoilc{ct},zzcoilc{ct},Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,i,10,Ncc);

trisurf(t2psym,psym(:,1),psym(:,2),psym(:,3),Att1,'facecolor','interp')


faces=[t2psym(:,1) t2psym(:,2);...
    t2psym(:,1) t2psym(:,3);...
    t2psym(:,2) t2psym(:,3);];
faces=sort(faces,2);
[foo,ix,jx]=unique(faces,'rows','stable');
vec=histc(jx,1:max(jx));
qx=find(vec==1);
ed=faces(ix(qx),:);
bpts=unique(ed(:));
close all
 trisurf(t2psym,psym(:,1),psym(:,2),psym(:,3),Att1,'edgealpha',0,'facecolor','interp')
 Att1=Att1/max(Att1);
 [t2p,tind]=TEMPgnerateholes(psym,t2psym,Att1,Fsym1,[.07 -.07])
 close all
 nrr=1:numel(t2p);
for i=nrr
rr(i)=numel(t2p{i});
end
[val,~,rr2]=unique(rr);

reg{1}=nrr(rr2==3);
reg{2}=nrr(rr2==2);
reg{3}=nrr(rr2==1);
for j=1:3
    t2p2{j}=t2p{reg{j}(1)};
    tind2{j}=tind{reg{j}(1)};
    for i=2:numel(reg{j})
        t2p2{j}=cat(1,t2p2{j},t2p{reg{j}(i)});
        tind2{j}=cat(1,tind2{j},tind{reg{j}(i)});

    end
end
%%
t2p=t2p2; 
for i=1:numel(t2p)
 %subplot(3,3,i),
hold on
col=i*ones([numel(t2p{i})/3 1]); 
trisurf(t2p{i},psym(:,1),psym(:,2),psym(:,3),col,'edgealpha',0,'facealpha',.4)
text(psym(t2p{i}(1,1),1),psym(t2p{i}(1,1),2),psym(t2p{i}(1,1),3),num2str(i)); 
end

load(harmfile)
load(designfile);
j=2;
Jx=Ax(:,1:Ncc)*Jcalc2{i};
Jy=Ay(:,1:Ncc)*Jcalc2{i};
Jz=Az(:,1:Ncc)*Jcalc2{i};
quiver3(loc(tind2{j},1),loc(tind2{j},2),loc(tind2{j},3),...
    Jx(tind2{j}),Jy(tind2{j}),Jz(tind2{j}));

