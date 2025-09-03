%% Step 1: generate and plot mesh
clear all
[x,y,z]=ndgrid(0:.016:.16, 0:.016:.16,.09);
tri=delaunay(x,y);
p(:,1)=x(:); p(:,2)=y(:); p(:,3)=z(:);
trisurf(tri,p(:,1),p(:,2),p(:,3));
view(2)
clear x y z
%% Step 2: generate mesh data structures
addpath('./generateharmonics')
harmfile='tutorial_harmonics';
nlayers=1;
debug=1;
nke=30; %number of modes to store
generateharmonicsinputmeshtri(harmfile,nlayers,p,tri,debug,nke);
%% step 3: Run optimization
addpath('./optimize');
Ncc=10;%number of harmonics should be less than nke
depth=140;%depth in .1mm
energy=[60 100 150 200 250];%different energy levels
filename='testoptim'; %coil design filename
filenamesav=strcat(filename,'d',num2str(depth),'ncoils',num2str(Ncc));
designfile=strcat(filenamesav,'lmax30.mat');
ncoils = runoptiminoga(depth,energy,Ncc,harmfile,filename);
%% see tradeoffs
addpath('.\determinetradeoff');
[Wi,Si,peakJ]=anycoildeterminetradeoff_f(harmfile,designfile,140,Ncc);
%% step 4: Discretize coil
%get windings
Nwindings=10;
designid=3;
nlayers=1;
[xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,nlayers,designid,Nwindings,Ncc);
close all%plot windings
%resample coil at .1mm per segment and save
for i=1:numel(xcoil)
    hold on
    if sign(zzcoil{i})==1
plot3(xcoil{i},ycoil{i},zcoil{i},'blue','linewidth',2);
    else
plot3(xcoil{i},ycoil{i},zcoil{i},'red','linewidth',2);
    end
    
    
 for i=1:numel(xcoil)
    len=sqrt((xcoil{i}(2:1:end)-xcoil{i}(1:end-1)).^2+...
             (ycoil{i}(2:1:end)-ycoil{i}(1:end-1)).^2+...
             (zcoil{i}(2:1:end)-zcoil{i}(1:end-1)).^2);
         len=cumsum([0;len(:)]);
         len(end)
         N=ceil(len(end)/.0001);
    locs=0:len(end)/N:len(end);
     xcoil{i}=interp1(len,xcoil{i},locs,'linear');
     ycoil{i}=interp1(len,ycoil{i},locs,'linear');
     zcoil{i}=interp1(len,zcoil{i},locs,'linear');    
end
save('..\tutorialcoil.mat','xcoil','ycoil','zcoil','zzcoil','Fsym1')
end
axis equal
view(2)
%% step 5 plan out windings
%this only works with coils that are topologically equivalent to the
%standard fdtms coil.
addpath('./generatecoilmesh');
filename='./newcoilsp.mat';
outfile='testouts';
plansinglewindingside2(filename,outfile);

load(strcat(outfile,'_fdtmscoilsp.mat'),'rv','jv','tri','p','Jdc','R','Ldc');

%%
axis equal
Ldc=real(Ldc);
addpath('./ANAT_Sphere_code/');
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
V=nnz(Ema>=50)*(.1)^3;
S=V/D;
W=(cur/(2*pi*3000))^2*real(Ldc)/2;
W
Ldc
R
S
D
Z=R+1i*2*pi*3000*Ldc