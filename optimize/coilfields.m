function [rs,js,Eout,p,t2p,clen,...
    xcoil,ycoil,zcoil,zzcoil,boo]=coilfields(p,t2p,Att,levs)


[cou,hout] = tricontour(p(:,1:2),t2p,Att,levs);
Fsym1=scatteredInterpolant(p(:,1),p(:,2),p(:,3),'natural','none');

close all
    ct=1;
    en=0;
c=cou;
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
zcoil{ct}=Fsym1(c(1,st:en),c(2,st:en));


zzcoil{ct}=height;
ct=ct+1;

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

clen=0;
boo=(numel(xcoil)>10*numel(levs));

for i=1:numel(xcoil)
    len=sqrt((xcoil{i}(2:1:end)-xcoil{i}(1:end-1)).^2+...
             (ycoil{i}(2:1:end)-ycoil{i}(1:end-1)).^2+...
             (zcoil{i}(2:1:end)-zcoil{i}(1:end-1)).^2);         
     len=cumsum([0;len(:)]);
     clen=clen+len(end);
     N=ceil(len(end)/.0008);
     locs=0:len(end)/N:len(end);
     xcoil{i}=interp1(len,xcoil{i},locs,'linear');
     ycoil{i}=interp1(len,ycoil{i},locs,'linear');
     zcoil{i}=interp1(len,zcoil{i},locs,'linear');    
end

if numel(xcoil)==numel(levs) || numel(xcoil)==2*numel(levs) || numel(xcoil)>=3*numel(levs)
else
nrem=numel(xcoil)-numel(levs);
if nrem>numel(levs)
nrem=numel(xcoil)-2*numel(levs);
end
for i=1:numel(xcoil)
lenii(i)=numel(xcoil{i});
end
[lenii,Ivv] = sort(lenii);
ct=1;
Ivv
for i=nrem+1:numel(xcoil)
xcoil2{ct}=xcoil{Ivv(i)};
ycoil2{ct}=ycoil{Ivv(i)};
zcoil2{ct}=zcoil{Ivv(i)};
ct=ct+1;
end

xcoil=xcoil2;ycoil=ycoil2;zcoil=zcoil2;
end
i=1;
clear rs js
rs(:,1)=cat(1,(xcoil{i}(2:end)+xcoil{i}(1:end-1))'/2,...
    (xcoil{i}(1)+xcoil{i}(end))/2);
rs(:,2)=cat(1,(ycoil{i}(2:end)+ycoil{i}(1:end-1))'/2,...
    (ycoil{i}(1)+ycoil{i}(end))/2);
rs(:,3)=cat(1,(zcoil{i}(2:end)+zcoil{i}(1:end-1))'/2,...
    (zcoil{i}(1)+zcoil{i}(end))/2);
js(:,1)=cat(1,(xcoil{i}(2:end)-xcoil{i}(1:end-1))',...
    (xcoil{i}(1)-xcoil{i}(end)));
js(:,2)=cat(1,(ycoil{i}(2:end)-ycoil{i}(1:end-1))',...
    (ycoil{i}(1)-ycoil{i}(end)));
js(:,3)=cat(1,(zcoil{i}(2:end)-zcoil{i}(1:end-1))',...
    (zcoil{i}(1)-zcoil{i}(end)));

for i=2:numel(xcoil)
    rst(:,1)=cat(1,rs(:,1),(xcoil{i}(2:end)+xcoil{i}(1:end-1))'/2,...
    (xcoil{i}(1)+xcoil{i}(end))/2);
    rst(:,2)=cat(1,rs(:,2),(ycoil{i}(2:end)+ycoil{i}(1:end-1))'/2,...
    (ycoil{i}(1)+ycoil{i}(end))/2);
    rst(:,3)=cat(1,rs(:,3),(zcoil{i}(2:end)+zcoil{i}(1:end-1))'/2,...
    (zcoil{i}(1)+zcoil{i}(end))/2);
    jst(:,1)=cat(1,js(:,1),(xcoil{i}(2:end)-xcoil{i}(1:end-1))',...
    (xcoil{i}(1)-xcoil{i}(end)));
    jst(:,2)=cat(1,js(:,2),(ycoil{i}(2:end)-ycoil{i}(1:end-1))',...
    (ycoil{i}(1)-ycoil{i}(end)));
    jst(:,3)=cat(1,js(:,3),(zcoil{i}(2:end)-zcoil{i}(1:end-1))',...
    (zcoil{i}(1)-zcoil{i}(end)));
rs=rst;js=jst;
clear rst jst
end
load barisamp4.mat
addpath('C:\Users\ljg24\Desktop\FEM_modes\matlab');
p=.073*p;
[Eout]=computeEprimary(rs',js',numel(rs)/3,p',numel(p)/3);


end


function iscw=ispolycw(x,y)
iscw=sign(sum((x(2:end)-x(1:end-1)).*(y(2:end)+y(1:end-1))));
iscw(iscw~=1)=0;

end