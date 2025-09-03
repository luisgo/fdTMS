clear all
harmfile='../newwireharmonics.mat';
dept=170;Wc=107;Sc=15.2;
Ncc=16;
coilid=4;
designfile=strcat('optimT3no100d',num2str(dept),'ncoils16lmax30.mat');
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    ct=1;
i=4;
%[xcoilbig,ycoilbig,zcoilbig,zzcoilbig,...
 %         xcoilfig8,ycoilfig8,zcoilfig8,zzcoilfig8,...
 %         xcoilside,ycoilside,zcoilside,zzcoilside]=generatecoilwindingsmultipolesallhighdepth(harmfile,designfile,1,coilid,8,Ncc);
 
 
[xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,coilid,10,Ncc);
close all
 %%
 for i=1:numel(xcoil)
     plot3(xcoil{i},ycoil{i},zcoil{i})
     hold on
 end
      %%
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

outfile='hh';
addpath('C:\Users\ljg24\Desktop\Coiloptcode\generatecoilmesh');
separateloops3(xcoilbig,ycoilbig,zcoilbig,...
          xcoilfig8,ycoilfig8,zcoilfig8,...
          xcoilside,ycoilside,zcoilside,outfile);
      
      plansinglewindingside(filename,outfile);     
      %
load(strcat(outfile,'_fdtmscoilsp.mat'),'rv','jv','tri','p','Jdc','R','Ldc');



%%
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
