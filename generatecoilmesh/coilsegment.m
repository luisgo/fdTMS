function [rs,nrs,stw,enw]=coilsegment(Vx,Vy,Vz,c1,c1ran,layer,layerp,rsp,res,Fsym,Vxt,Vyt,Vzt,nrsp);
load firstordercorr_middle.mat nx ny nz
cent(1)=mean(Vxt{c1(1)});
cent(2)=mean(Vyt{c1(1)});
cent(3)=mean(Vzt{c1(1)});
cen(1)=mean(Vx{c1(1)});
cen(2)=mean(Vy{c1(1)});
cen(3)=mean(Vz{c1(1)});
cen=cen-cent;
len=cumsum(sqrt((Vx{c1(1)}(2:end)-Vx{c1(1)}(1:end-1)).^2+...
                (Vy{c1(1)}(2:end)-Vy{c1(1)}(1:end-1)).^2+...
                (Vz{c1(1)}(2:end)-Vz{c1(1)}(1:end-1)).^2));
len=[0;len(:)];
rang=min(len):(max(len)-min(len))/numel(Vxt{c1(1)}(:)):max(len);
[~,uind,~]=unique(len);

Vx{c1(1)}=interp1(len(uind),Vx{c1(1)}(uind),rang,'linear');
Vy{c1(1)}=interp1(len(uind),Vy{c1(1)}(uind),rang,'linear');
Vz{c1(1)}=interp1(len(uind),Vz{c1(1)}(uind),rang,'linear');
[~,I]=sort(sqrt((Vxt{c1(1)}(c1ran(1))-Vx{c1(1)}(:)+cen(1)).^2+...
                (Vyt{c1(1)}(c1ran(1))-Vy{c1(1)}(:)+cen(2)).^2+...
                (Vzt{c1(1)}(c1ran(1))-Vz{c1(1)}(:)+cen(3)).^2),1);
c1ran(1)
I(1)
            beg=I(1);
 
[~,I]=sort(sqrt((Vxt{c1(1)}(c1ran(end))-Vx{c1(1)}(:)+cen(1)).^2+...
                (Vyt{c1(1)}(c1ran(end))-Vy{c1(1)}(:)+cen(2)).^2+...
                (Vzt{c1(1)}(c1ran(end))-Vz{c1(1)}(:)+cen(3)).^2),1);
            
ending=I(1);
if ending<beg
c1ran=[beg:numel(Vx{c1(1)}) 1:ending];
else
c1ran=beg:ending;
end

load layerinfo.mat layerp;
rs(:,1)=Vx{c1(1)}(c1ran(:))+(layer)*nx(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,2)=Vy{c1(1)}(c1ran(:))+(layer)*ny(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,3)=Vz{c1(1)}(c1ran(:))+(layer)*nz(Vx{c1(1)}(c1ran(:)),Vy{c1(1)}(c1ran(:)));
rs(:,4)=Vx{c1(1)}(c1ran(:));
rs(:,5)=Vy{c1(1)}(c1ran(:));

nrs=floor(numel(rs)/8);

stw=1;
enw=numel(rs(:,1));
if numel(rsp)~=0

len=sqrt(sum((rsp(end,:)-rs(1,:)).^2,2));
N=ceil(len/res);
co=1/N:1/N:1-1/N;

rstemp(:,1)=(rsp(end,1)-layerp*nx(rsp(end,1),rsp(end,2)))*(1-co(:))+Vx{c1(1)}(c1ran(1))*(co(:));
rstemp(:,2)=(rsp(end,2)-layerp*ny(rsp(end,1),rsp(end,2)))*(1-co(:))+Vy{c1(1)}(c1ran(1))*(co(:));
rstemp(:,3)=Fsym(rstemp(:,1),rstemp(:,2));
rst(:,1)=rstemp(:,1)+(layerp*(1-co(:))+layer*co(:)).*nx(rstemp(:,1),rstemp(:,2));
rst(:,2)=rstemp(:,2)+(layerp*(1-co(:))+layer*co(:)).*ny(rstemp(:,1),rstemp(:,2));
rst(:,3)=rstemp(:,3)+(layerp*(1-co(:))+layer*co(:)).*nz(rstemp(:,1),rstemp(:,2));
rst(:,4)=rstemp(:,1);
rst(:,5)=rstemp(:,2);
Npts=numel(rst(:,1));
stid=numel(rsp(:,1));

 Npad=150;
if Npad>nrs || Npad>nrsp
     Npad=min(nrs,nrsp);
end
ending=stid+Npts+Npad;
rs=cat(1,rsp,rst,rs);
stw=numel(rsp(:,1))+numel(rst(:,1));
enw=numel(rs(:,1));
if stid-Npad>0
     rs(stid-Npad:ending,1)= movmean(rs(stid-Npad:ending,1),Npad);
 rs(stid-Npad:ending,2)= movmean(rs(stid-Npad:ending,2),Npad);
 rs(stid-Npad:ending,3)= movmean(rs(stid-Npad:ending,3),Npad);
 rs(stid-Npad:ending,4)= movmean(rs(stid-Npad:ending,4),Npad);
rs(stid-Npad:ending,5)= movmean(rs(stid-Npad:ending,5),Npad);
end
end
layerp=layer;
save layerinfo.mat layerp;
end