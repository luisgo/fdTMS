function [t2p,tind,Attmod,pmod,bdr]=TEMPgnerateholesnew(psym,t2psym,Att1,Fsym1,levs)
 [xcoil,ycoil,zcoil,zzcoil]=getcontours(psym,t2psym,Fsym1,Att1,levs);
close all
 [xzero,yzero,zzero,zzzero]=getcontours(psym,t2psym,Fsym1,Att1,[abs(levs(1))-.015 -abs(levs(1))+.015]);
 

 %reorient contours
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



TR = triangulation(t2psym, psym(:,1), psym(:,2), psym(:,3));
for ival=1:numel(xcoil)
close all
%get interior points to the poligon
[in,on] = inpolygon(psym(:,1),psym(:,2),xcoil{ival},ycoil{ival});
in(on==1)=0;
hold on
%get triangles interior to polygon
id=1:numel(in);
TI = vertexAttachments(TR, id(in)');
ct=0
for i=1:numel(TI)
ct=ct+numel(TI{i});
end
trian=zeros([ct 1]);
ct=0;
for i=1:numel(TI)
trian(ct+1:ct+numel(TI{i}))=TI{i};
    ct=ct+numel(TI{i});
end

t2p{ival}=t2psym(unique(trian),:);
tind{ival}=unique(trian);
bdrnode=setdiff(unique(t2p{ival}(:)),id(in));
pmod{ival}=psym;
Attmod{ival}=Att1;
scorr=sign(mean(Att1(id(in))));
DT=triangulation(t2p{ival},pmod{ival}(:,1),pmod{ival}(:,2));
for i=1:numel(bdrnode)
    dist=10000;
    if psym(bdrnode(i),2)==0
        pmod{ival}(bdrnode(i),:)=psym(bdrnode(i),:);
    
               pmod{ival}(bdrnode(i),2)=sign(scorr)*0;
        pmod{ival}(bdrnode(i),3)=Fsym1(pmod{ival}(bdrnode(i),1),pmod{ival}(bdrnode(i),2));
        [TI, BC] = pointLocation(DT, pmod{ival}(bdrnode(i),1:2));
        ou=Attmod{ival}(t2p{ival}(TI,1))*BC(1)+...
           Attmod{ival}(t2p{ival}(TI,2))*BC(2)+...
           Attmod{ival}(t2p{ival}(TI,3))*BC(3)
       
       Attmod{ival}(bdrnode(i))= ou;
       continue;
    end
for j=1:numel(xzero)
    if sign(zzzero{j})~=scorr
        continue;
    end
    for k=1:numel(xzero{j})-1
    tan=[xzero{j}(k+1)-xzero{j}(k) yzero{j}(k+1)-yzero{j}(k)];
    
    u=[xzero{j}(k)-psym(bdrnode(i),1) yzero{j}(k)-psym(bdrnode(i),2)];
    t=-sum(u.*tan)/sum(tan.*tan);
    if t>=0 && t<=1
        ou2(1)=(1-t)*xzero{j}(k)+t*xzero{j}(k+1);
        ou2(2)=(1-t)*yzero{j}(k)+t*yzero{j}(k+1);
        no2=norm(psym(bdrnode(i),1:2)'-ou2(:));
    else
        no2=norm(psym(bdrnode(i),1:2)-[xzero{j}(k) yzero{j}(k)]);
            ou2=[xzero{j}(k) yzero{j}(k)];
        no3=norm(psym(bdrnode(i),1:2)-[xzero{j}(k+1) yzero{j}(k+1)]);
        if no3<no2
            ou2=[xzero{j}(k+1) yzero{j}(k+1)];
            no2=no3;
        end
    end
        if dist>no2 
           pmod{ival}(bdrnode(i),1:2)=ou2;
           pmod{ival}(bdrnode(i),3)=Fsym1(ou2(1),ou2(2));
           Attmod{ival}(bdrnode(i))=zzzero{j};
           dist=no2;
        end

    
    end
end
end
t2ptemp=cat(1,t2p{ival}(:),bdrnode(:));
[pp,~,t2ptemp]=unique(t2ptemp);
t2p{ival}=reshape(t2ptemp(1:numel(t2p{ival})),size(t2p{ival}));
pmod{ival}=pmod{ival}(pp,:);
Attmod{ival}=Attmod{ival}(pp);
bdr{ival}=t2ptemp(numel(t2p{ival})+1:end);



end
end


function [xcoil,ycoil,zcoil,zzcoil]=getcontours(psym,t2psym,Fsym1,Att1,levs)
 [cou{1},hout] = tricontour(psym(:,1:2),t2psym,Att1,levs); 
 
      ct=1;
for i=1:1
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

end
end

end