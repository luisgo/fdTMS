function [t2p,tind]=TEMPgnerateholes(psym,t2psym,Att1,Fsym1,levs)
 [xcoil,ycoil,zcoil,zzcoil]=getcontours(psym,t2psym,Fsym1,Att1,levs);
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
% dd=t2psym(unique(trian),:);
% dd=unique(dd(:));
% %get triangles interior to polygon
% in(dd)=1;
% id=1:numel(in);
% TI = vertexAttachments(TR, id(in)');
% ct=0
% for i=1:numel(TI)
% ct=ct+numel(TI{i});
% end
% trian=zeros([ct 1]);
% ct=0;
% for i=1:numel(TI)
% trian(ct+1:ct+numel(TI{i}))=TI{i};
%     ct=ct+numel(TI{i});
% end

t2p{ival}=t2psym(unique(trian),:);
tind{ival}=unique(trian);
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