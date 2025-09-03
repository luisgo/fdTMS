function [xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,nlayers,coilid,N,Ncc)
load(harmfile);
load(designfile);

psym2=psym;
[pind,~,tp]=unique(t2psym(:));
psym=psym(pind,:);
norm(psym2(:)-psym(:))
t2psym=reshape(tp,size(t2psym));


if nlayers==2
faces=[t2psym(:,[1,2]);
       t2psym(:,[1,3]);
       t2psym(:,[2,3])];
[foo,ix,jx]=unique(faces,'rows','stable');
faces=faces(ix,:);
G = graph(faces(:,1),faces(:,2));
BINS = conncomp(G);
nod=1:numel(psym(:,1));
ran1p=nod(BINS==1);
ran2p=nod(BINS==2);
Fsym1=scatteredInterpolant(psym(ran1p,1),psym(ran1p,2),psym(ran1p,3),'natural','none');
Fsym2=scatteredInterpolant(psym(ran2p,1),psym(ran2p,2),psym(ran2p,3),'natural','none');
nt=numel(t2psym(:,1))/4;
ran1=cat(1,1:1*nt/2,nt+1:3*nt/2,2*nt+1:5*nt/2,3*nt+1:7*nt/2);
ran2=cat(1,round(nt/2)+1:nt,3*round(nt/2)+1:2*nt,5*round(nt/2)+1:3*nt,7*round(nt/2)+1:4*nt);

else
Fsym1=scatteredInterpolant(psym(:,1),psym(:,2),psym(:,3),'natural','none');
Fsym2=[];
end
%scatter3(psym(ran1p,1),psym(ran1p,2),psym(ran1p,3))
%hold on
%scatter3(psym(ran2p,1),psym(ran2p,2),psym(ran2p,3))

Att=Asym*Atrans(:,1:Ncc)*Jcalc2{coilid};
Att=Att(pind,:);


%%
clear xcoil ycoil zcoil


ms=max(Att);
mi=min(Att);
del=(ms-mi)/N;
levs=mi+del/2:del:ms-del/2;

Att1=Att;Att2=Att;
if nlayers==2
Att1(ran2p)=0;
Att2(ran1p)=0;
[cou{1},hout] = tricontour(psym(:,1:2),t2psym(ran1,:),Att1,levs);
[cou{2},hout] = tricontour(psym(:,1:2),t2psym(ran2,:),Att2,levs);
else
[cou{1},hout] = tricontour(psym(:,1:2),t2psym,Att1,levs);
end
save coilforplay psym t2psym Att1 Fsym1;

close all

%subplot(1,3,1),
view([0 0 1]);
hold on


    ct=1;
for i=1:nlayers
    if numel(cou{i})==0
        continue;
    end
    en=0;
c=cou{i};
while en<numel(c(1,:))
    st=en+1;
en=c(2,st)+st;
height=c(1,st);
if sign(c(1,st))==1
    col='r';
    sig(ct)=1;
else
    sig(ct)=-1;
    col='b';
end
st=st+1;

clear temp

xcoil{ct}=c(1,st:en);
ycoil{ct}=c(2,st:en);
if i==1
zcoil{ct}=Fsym1(c(1,st:en),c(2,st:en));
else
zcoil{ct}=Fsym2(c(1,st:en),c(2,st:en));    
end
zzcoil{ct}=height;
ct=ct+1;
if i==1
plot3(c(1,st:en),c(2,st:en),Fsym1(c(1,st:en),c(2,st:en)),'blue','linewidth',2);
elseif i==2
plot3(c(1,st:en),c(2,st:en),Fsym2(c(1,st:en),c(2,st:en)),'red','linewidth',2);
end
%scatter3(c(1,st:en),c(2,st:en),Fsym(c(1,st:en),c(2,st:en)))
hold on

end
end
for i=1:numel(xcoil);
if ispolycw(xcoil{i}, ycoil{i})==1 && sign(zzcoil{i})==1
elseif ispolycw(xcoil{i}, ycoil{i})==0 && sign(zzcoil{i})==-1
else
    xcoil{i}=xcoil{i}(end:-1:1);
    ycoil{i}=ycoil{i}(end:-1:1);
    zcoil{i}=zcoil{i}(end:-1:1);
end
end
axis equal
end
function iscw=ispolycw(x,y)
iscw=sign(sum((x(2:end)-x(1:end-1)).*(y(2:end)+y(1:end-1))));
iscw(iscw~=1)=0;

end