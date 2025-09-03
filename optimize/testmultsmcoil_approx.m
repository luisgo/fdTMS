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
wid=4.067098615315901e-03/2;
[xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,coilid,6,Ncc);
load firstordercorr_middle.mat nx ny nz
rs=[]
js=[]
for i=1:numel(xcoil)
    hold on
    clear nhat
    nhat(:,1)=nx(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
    nhat(:,2)=ny(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
    nhat(:,3)=nz(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
          plot3(xcoil{i}([1:end ])+wid*nhat(:,1),...
                ycoil{i}([1:end ])+wid*nhat(:,2),...
                zcoil{i}([1:end ])+wid*nhat(:,3),'r','linewidth',2)
          plot3(xcoil{i}([1:end ])-wid*nhat(:,1),...
                ycoil{i}([1:end ])-wid*nhat(:,2),...
                zcoil{i}([1:end ])-wid*nhat(:,3),'b','linewidth',2)
 

ra=[xcoil{i}([1:end ])'+wid*nhat(:,1),...
ycoil{i}([1:end ])'+wid*nhat(:,2),...
zcoil{i}([1:end ])'+wid*nhat(:,3)];
rb=[xcoil{i}([1:end ])'-wid*nhat(:,1),...
ycoil{i}([1:end ])'-wid*nhat(:,2),...
zcoil{i}([1:end ])'-wid*nhat(:,3)];      
len=cumsum(sqrt([0;sum((ra(2:end,:)-ra(1:end-1,:)).^2,2)]));
xa = interp1(len,ra(:,1),0:.0002:len(end));
ya = interp1(len,ra(:,2),0:.0002:len(end));
za = interp1(len,ra(:,3),0:.0002:len(end));
ra=[xa(:),ya(:),za(:)];

len=cumsum(sqrt([0;sum((rb(2:end,:)-rb(1:end-1,:)).^2,2)]));
xb = interp1(len,rb(:,1),0:.0002:len(end));
yb = interp1(len,rb(:,2),0:.0002:len(end));
zb = interp1(len,rb(:,3),0:.0002:len(end));
rb=[xb(:),yb(:),zb(:)];


rs=cat(1,rs,...
    (ra(2:end,:)+ra(1:end-1,:))/2,...
    (rb(2:end,:)+rb(1:end-1,:))/2);       
js=cat(1,js,...
    (ra(2:end,:)-ra(1:end-1,:)),...
    (rb(2:end,:)-rb(1:end-1,:)));
end
close all
quiver3(rs(:,1),rs(:,2),rs(:,3),...
    js(:,1),js(:,2),js(:,3))

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
