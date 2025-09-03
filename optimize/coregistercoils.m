path('C:\Users\ljg24\Desktop\Coiloptcode\generatecoilmesh\figure8coilmodelgenerator\wirecode',path);
path('C:\Users\ljg24\Desktop\Coiloptcode\generatecoilmesh\figure8coilmodelgenerator',path);
%%
x=load('C:\Users\ljg24\Desktop\Coiloptcode\Magstim_70mm_Fig8.ccd');
Ef=niftiread('C:\Users\ljg24\Desktop\Coiloptcode\Magstim_70mm_Fig8.nii.gz');
info=niftiinfo('C:\Users\ljg24\Desktop\Coiloptcode\Magstim_70mm_Fig8.nii.gz');
addpath C:\Users\ljg24\Desktop\FEM_modes\matlab
%%

[X Y Z]=ndgrid(1:121,1:81,1:41);
X=.005*(X-61);
Y=.005*(Y-41);
Z=.005*(Z-1);
robs=[X(:) Y(:) Z(:)]';
[Hout]=computeHprimary(x(:,1:3)',x(:,4:6)',numel(x)/6,robs,numel(robs)/3);
n=size(Ef);
Hout=reshape(real(Hout)',n);


%%
[Jdc,rv,Ldc,R,te2p,p]=figure8coilmodel1layerparallel(31,.023,.046,.015,.0003,.046,.002,1);
rv(:,3)=-rv(:,3)-.0035;
Jdc(:,3)=-Jdc(:,3);
p(:,3)=-p(:,3)-.0035;
t2p=surftri(p,te2p);
trisurf(t2p,p(:,1),p(:,2),p(:,3))

n=size(Ef);
[Eout]=computeEprimary(rv',real(Jdc'),numel(rv)/3,robs,numel(robs)/3);
Eout=reshape(real(Eout)',n);


niftiwrite(single(-Ef2), 'B35.nii',info);
%%
Ef2=niftiread('B35.nii');
subplot(1,2,1),
imagesc(Ef(:,:,20,2))

subplot(1,2,2),
imagesc(Ef2(:,:,20,2))


%%
[Jdc,rv,Ldc,R,te2p,p]=figure8coilmodel2layer(5,.035,.075,.006,.003,.078,.002,1);
rv(:,3)=-rv(:,3)-.0035;
Jdc(:,3)=-Jdc(:,3);
p(:,3)=-p(:,3)-.0035;
t2p=surftri(p,te2p);
trisurf(t2p,p(:,1),p(:,2),p(:,3))

n=size(Ef);
[Eout]=computeEprimary(rv',real(Jdc'),numel(rv)/3,robs,numel(robs)/3);
Eout=reshape(real(Eout)',n);

niftiwrite(single(Eout), 'B65.nii',info);
%%
Ef2=niftiread('B65.nii');
subplot(1,2,1),
imagesc(Ef(:,:,20,2))

subplot(1,2,2),
imagesc(Ef2(:,:,20,2))


%%
[Jdc,rv,Ldc,~,te2p,p]=figure8coilmodel2layerdouble([4 3],.067,(.095-.067)/8,.0059,.003,.100,.002,1);
rv(:,3)=-rv(:,3)-.002;
Jdc(:,3)=-Jdc(:,3);
p(:,3)=-p(:,3)-.002;
t2p=surftri(p,te2p);
trisurf(t2p,p(:,1),p(:,2),p(:,3))

n=size(Ef);
[Eout]=computeEprimary(rv',real(Jdc'),numel(rv)/3,robs,numel(robs)/3);
Eout=reshape(real(Eout)',n);

niftiwrite(single(Eout), 'B80.nii',info);
%%
Ef2=niftiread('B80.nii');
subplot(1,2,1),
imagesc(squeeze(Ef(:,41,:,2)))

subplot(1,2,2),
imagesc(squeeze(Ef2(:,41,:,2)))

%%
load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\midcoil_fdtmscoilsp.mat'),'rv','jv','tri','p','Jdc','R','Ldc');

rv=rv(:,[2,1,3]);
Jdc=Jdc(:,[2,1,3]);
p=p(:,[2,1,3]);
rv(:,3)=-(rv(:,3)-.087);
Jdc(:,3)=-Jdc(:,3);
p(:,3)=-(p(:,3)-.087);



n=size(Ef);
[Eout]=computeEprimary(rv',real(Jdc'),numel(rv)/3,robs,numel(robs)/3);


Eout=reshape(real(Eout)',n);


niftiwrite(single(-Eout), 'fdtmsb65.nii',info);
Ef2=niftiread('fdtmsb65.nii');

subplot(1,2,1),
imagesc(squeeze(Ef(:,:,10,2)))

subplot(1,2,2),
imagesc(squeeze(Ef2(:,:,10,2)))



%%
load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\midcoil_fdtmscoilsp.mat'),'rv','jv','tri','p','Jdc','R','Ldc');
rv=rv(:,[2,1,3]);
Jdc=Jdc(:,[2,1,3]);
p=p(:,[2,1,3]);
rv(:,3)=-(rv(:,3)-.087);
Jdc(:,3)=-Jdc(:,3);
p(:,3)=-(p(:,3)-.087);

addpath('../../reciprocitycodes/femcodes');
load ../../coiloptcodesphereres.mat p te2p;
te2p=te2p';p=p';
t2p=surftri(p,te2p);
trisurf(t2p,p(:,1),p(:,2),p(:,3))

p(:,3)=p(:,3)-.087;
rv(:,1)=-rv(:,1);
rv(:,3)=-rv(:,3)+.087;
Jdc(:,1)=-Jdc(:,1);
Jdc(:,3)=-Jdc(:,3);
conductivity=ones([numel(te2p)/4 1]);
[solution]=currfem_cl(te2p',p',conductivity,rv,real(Jdc),1,1,1);
load '../../../Downloads/sphere_homogen_notrefined_PointsAreInGM.mat'
Eout2=FEM_interpolator(solution,Field1.node'/1000);
addpath('C:\Users\ljg24\Desktop\Coiloptcode\generatecoilmesh\figure8coilmodelgenerator\ANAT_Sphere_code')
Efield=computefields(rv,real(Jdc),zeros(size(Jdc)),zeros([numel(Jdc)/3 1]),1,25,.086,Field1.node(:,1:10^6)'/1000);