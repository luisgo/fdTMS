clear all;
harmfile='../newwireharmonics.mat';
addpath C:\Users\ljg24\Desktop\Coiloptcode\generatecoilmesh\;

dept=101;Wc=107;Sc=15.2;
Ncc=24;
coilid=6;
designfile=strcat('optimT3no100d',num2str(dept),'ncoils24lmax30.mat');
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Si','Wi');
ct=1;
i=5;
[xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,nlayers,desid,10,Ncc);
load(filename,'xcoil','ycoil','zcoil','Fsym1');
plansinglewinding(designfile,'trysmall.mat')
%[xcoilbig, ycoilbig, zcoilbig, zzcoilbig,...
% xcoilfig8,ycoilfig8,zcoilfig8,zzcoilfig8,...
% xcoilside,ycoilside,zcoilside,zzcoilside]=generatecoilwindingsmultipolesall2(harmfile,designfile,1,coilid,10,Ncc);
%close all
%save smallcoil3.mat
 %%
 load smallcoil2.mat
 for i=1:numel(xcoilfig8)
     plot3(xcoilfig8{i},ycoilfig8{i},zcoilfig8{i})
     hold on
       hold on
     text(xcoilfig8{i}(end),ycoilfig8{i}(end),zcoilfig8{i}(end),num2str(i))
 end
 for i=1:numel(xcoilbig)
     plot3(xcoilbig{i},ycoilbig{i},zcoilbig{i})
     hold on
     text(xcoilbig{i}(end),ycoilbig{i}(end),zcoilbig{i}(end),num2str(i))
 end
 %
  for i=1:numel(xcoilside)
     plot3(xcoilside{i},ycoilside{i},zcoilside{i})
     hold on
       hold on
     text(xcoilside{i}(end),ycoilside{i}(end),zcoilside{i}(end),num2str(i))
 end
      %
      ct=1;
      clear xcoil ycoil zcoil
      for i=1:numel(xcoilbig)
xcoil{ct}=xcoilbig{i};
ycoil{ct}=ycoilbig{i};
zcoil{ct}=zcoilbig{i};
ct=ct+1;
      end
clen=0;
resolu=.0001;
for i=1:numel(xcoil)
    len=sqrt((xcoil{i}(2:1:end)-xcoil{i}(1:end-1)).^2+...
             (ycoil{i}(2:1:end)-ycoil{i}(1:end-1)).^2+...
             (zcoil{i}(2:1:end)-zcoil{i}(1:end-1)).^2);    
         
     len=cumsum([0;len(:)]);
     clen=clen+len(end);
     N=ceil(len(end)/resolu);
     locs=0:len(end)/N:len(end);
     xcoil{i}=interp1(len,xcoil{i},locs,'linear');
     ycoil{i}=interp1(len,ycoil{i},locs,'linear');
     zcoil{i}=interp1(len,zcoil{i},locs,'linear');    
end
xcoilbig=xcoil;
ycoilbig=ycoil;
zcoilbig=zcoil;

xcoil=xcoilfig8;
ycoil=ycoilfig8;
zcoil=zcoilfig8;
clen=0;
for i=1:numel(xcoil)
    len=sqrt((xcoil{i}(2:1:end)-xcoil{i}(1:end-1)).^2+...
             (ycoil{i}(2:1:end)-ycoil{i}(1:end-1)).^2+...
             (zcoil{i}(2:1:end)-zcoil{i}(1:end-1)).^2);         
     len=cumsum([0;len(:)]);
     clen=clen+len(end);
     N=ceil(len(end)/resolu);
     locs=0:len(end)/N:len(end);
     xcoil{i}=interp1(len,xcoil{i},locs,'linear');
     ycoil{i}=interp1(len,ycoil{i},locs,'linear');
     zcoil{i}=interp1(len,zcoil{i},locs,'linear');    
end
xcoilfig8=xcoil;
ycoilfig8=ycoil;
zcoilfig8=zcoil;



xcoil=xcoilside;
ycoil=ycoilside;
zcoil=zcoilside;
clen=0;
ct=1;
for i=1:numel(xcoil)
    len=sqrt((xcoil{i}(2:1:end)-xcoil{i}(1:end-1)).^2+...
             (ycoil{i}(2:1:end)-ycoil{i}(1:end-1)).^2+...
             (zcoil{i}(2:1:end)-zcoil{i}(1:end-1)).^2);         
     len=cumsum([0;len(:)]);
     clen=clen+len(end);
     N=ceil(len(end)/resolu);
     locs=0:len(end)/N:len(end);
     xcoili{ct}=interp1(len,xcoil{i},locs,'linear');
     ycoili{ct}=interp1(len,ycoil{i},locs,'linear');
     zcoili{ct}=interp1(len,zcoil{i},locs,'linear');    
     ct=ct+1;
end
xcoilside=xcoili;
ycoilside=ycoili;
zcoilside=zcoili;
%
outfile='hh';
addpath('C:\Users\ljg24\Desktop\Coiloptcode\generatecoilmesh');
% separateloops3(xcoilbig,ycoilbig,zcoilbig,...
%           xcoilfig8,ycoilfig8,zcoilfig8,...
%           xcoilside,ycoilside,zcoilside,outfile);
separateloopstwolayer(xcoilbig,ycoilbig,zcoilbig,...
          xcoilfig8,ycoilfig8,zcoilfig8,...
          xcoilside,ycoilside,zcoilside,outfile);
       
      %%
load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\'...
    ,outfile,'_fdtmscoilsp.mat'),'rv','jv','tri','p','Jdc','R','Ldc');



%
axis equal
Ldc=real(Ldc);
addpath('C:\Users\ljg24\Desktop\Coiloptcode\generatecoilmesh\figure8coilmodelgenerator\ANAT_Sphere_code');
[Xeval,Yeval,Zeval]=ndgrid(-.07:.001:.07,-.07:.001:.07,.05:.001:.07);
pt(:,1)=Xeval(:);
pt(:,2)=Yeval(:);
pt(:,3)=Zeval(:);
Efield=computefields(rv,Jdc,zeros(size(rv)),zeros(size(rv(:,1))),...
    1,25,.086,pt);

%
pt2(:,3)=.07:-.0001:.035;
pt2(:,2)=0;
pt2(:,1)=0;
Efield2=computefields(rv,Jdc,zeros(size(rv)),zeros(size(rv(:,1))),...
    1,25,.086,pt2);

%
Efield=real(Efield);
Efield2=real(Efield2);
plot(Efield2(:,2))
Ema2=sqrt(sum(Efield2.^2,2));
Ema=sqrt(sum(Efield.^2,2));
cur=100/max(Ema2);
Ema2=cur*Ema2;
Ema=cur*Ema;
r=sqrt(sum(pt.^2,2));
Ema(r>.07)=0;
D=(nnz(Ema2>=50)-1)/100;
if D*100>dept
Ema2=sqrt(sum(Efield2.^2,2));
Ema=sqrt(sum(Efield.^2,2));
cur=50/Ema2(dept+1);
Ema2=cur*Ema2;
Ema=cur*Ema;
r=sqrt(sum(pt.^2,2));
Ema(r>.07)=0;
D=(nnz(Ema2>=50)-1)/100;
 
end

V=nnz(Ema>=50)*(.1)^3;
S=V/D;
W=(cur/(2*pi*3000))^2*real(Ldc)/2;
W
Ldc
R
S
D
Z=R+1i*2*pi*3000*Ldc
%%
close all
quiver3(rv(:,1),rv(:,2),rv(:,3),Jdc(:,1),Jdc(:,2),Jdc(:,3))
axis equal
view(2)
