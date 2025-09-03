function [xcoil,ycoil,zcoil,zzcoil]=getcontours(psym,t2psym,Att1,Fsym1,levs)
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

for ival=1:numel(xcoil)

   len=sqrt((xcoil{ival}(2:1:end)-xcoil{ival}(1:end-1)).^2+...
            (ycoil{ival}(2:1:end)-ycoil{ival}(1:end-1)).^2+...
            (zcoil{ival}(2:1:end)-zcoil{ival}(1:end-1)).^2);
         len=cumsum([0;len(:)]);
         len(end)
         N=ceil(len(end)/.0001);
     
    locs=0:len(end)/N:len(end);
     xcoil{ival}=interp1(len,xcoil{ival},locs,'linear');
     ycoil{ival}=interp1(len,ycoil{ival},locs,'linear');
     zcoil{ival}=interp1(len,zcoil{ival},locs,'linear');    
end


 %reorient contours
for i=1:numel(xcoil);
if ispolycw(xcoil{i}, ycoil{i})==1 && sign(zzcoil{i})==1
elseif ispolycw(xcoil{i}, ycoil{i})==0 && sign(zzcoil{i})==-1
else
    xcoil{i}=xcoil{i}(end:-1:1);
    ycoil{i}=ycoil{i}(end:-1:1);
    zcoil{i}=zcoil{i}(end:-1:1);
end

hold on
text(xcoil{i}(1),ycoil{i}(1),zcoil{i}(1),num2str(i));
end
end