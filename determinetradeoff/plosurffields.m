load baricooroNEW4.mat
load 'newharmany.mat'

Xeval=Xeval(Zeval>0);
Yeval=Yeval(Zeval>0);
Zeval=Zeval(Zeval>0);
Xeval=cat(1,Xeval(:),-Xeval(:),Xeval(:),-Xeval(:));
Yeval=cat(1,Yeval(:),Yeval(:),-Yeval(:),-Yeval(:));
Zeval=cat(1,Zeval(:),Zeval(:),Zeval(:),Zeval(:));
tri=delaunay(Xeval,Yeval);
trisurf(tri,Xeval,Yeval,Zeval);
[Emat]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);

short=1;
%%
ncoils2=[15];
ndpt=[140];

for cc=1
fold='C:\Users\ljg24\Desktop\optimisationcodes\anycoils\';
rrri=strcat(fold,'fourdesigns',num2str(ndpt(cc)),'front20ncoils',num2str(ncoils2(cc)),'.mat');    

Ndep=ndpt(cc);
    clear W V
load(rrri)
ncoils=ncoils2(cc);
Lmat=Lmat2(1:ncoils,1:ncoils);
for des=1:numel(Jcalc2)
    if numel(Jcalc2{des})==0 || nnz(Jcalc2{des})==0
        continue;
    end
Jcalc=Jcalc2{des}(1:ncoils);
N=size(Jcalc);
Nharm=numel(Emat);
%populate field matrices
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
stack=globind(lmax,lmax);
coeff=zeros(size(harmcoeff{1,1}(:)));
for icoil=1:ncoils
coeff(:)=coeff(:)+Jcalc(icoil)*real(harmcoeff{icoil}(:)); %harmonic currents to spherical currents
end


E{des}=zeros(size(Emat{1}));
for l=1:lmax
for m=1:1:l
E{des}=E{des}+coeff(globind(m,l)+stack)*Emat{globind(m,l)+stack};
end
end

end
end
%%

load C:\Users\ljg24\Desktop\optimisationcodes\anycoils\surfpts.mat Xeval Yeval Zeval E tri;
load(strcat('mmS','coil.mat'));
load C:\Users\ljg24\Desktop\optimisationcodes\anycoils\surfpts.mat Xeval Yeval Zeval E tri;

rv(:,3)=rv(:,3)++.090;
n=size(Xeval);
pt(:,1)=Xeval(:);
pt(:,2)=Yeval(:);
pt(:,3)=Zeval(:);
    Efield=computefields(rv,Jdc,zeros(size(rv)),zeros(size(rv(:,1)))...
       ,2*pi*3000,30,.084,pt);
   con =3.7701e+03;

Ema=con*sqrt(sum(abs(Efield).^2,2));

 
save figure8fields.mat
%%
subplot(2,2,1),
trisurf(tri,Xeval,Yeval,Zeval,Ema,'edgealpha',0,'facecolor','interp')
view([0 0 1]);
axis equal
axis off



ax = gca;
mymap = colormap(ax);
newmap=zeros([128 3]);
newmap(:,1)=interp1(sqrt(1:64)',mymap(:,1),1:7/127:8);
newmap(:,2)=interp1(sqrt(1:64)',mymap(:,2),1:7/127:8);
newmap(:,3)=interp1(sqrt(1:64)',mymap(:,3),1:7/127:8);

colormap(newmap);


%%
pli=[7 5 3 1]
for i=1:4
Emag=sqrt(sum(abs(E{pli(i)}).^2,4));
subplot(2,2,i),
trisurf(tri,Yeval,Xeval,Zeval,(Emag)/max(Emag),'edgealpha',0,'facecolor','interp')
view([0 0 1]);
axis equal
axis off



ax = gca;
mymap = colormap(ax);
newmap=zeros([128 3]);
newmap(:,1)=interp1(sqrt(1:64)',mymap(:,1),1:7/127:8);
newmap(:,2)=interp1(sqrt(1:64)',mymap(:,2),1:7/127:8);
newmap(:,3)=interp1(sqrt(1:64)',mymap(:,3),1:7/127:8);

colormap(newmap);
end