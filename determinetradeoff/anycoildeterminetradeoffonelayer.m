%clear all
lmax=30;
short=1;
load 'newharmany.mat'
if short==1
    
load shortcut.mat Emat2 Emat;
else
[Yeval,Xeval,Zeval]=meshgrid(-.07:.001:.07,-.07:.001:.07,0.0305:.001:.07);
[Emat]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
[Yeval,Xeval,Zeval]=meshgrid(0,0,.07:-.0001:.01);
[Emat2]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
save shortcut.mat -v7.3 Emat2 Emat;
end
%%

[Yeval,Xeval,Zeval]=meshgrid(0,0,.07:-.0001:.01);
[Emat2]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
ncoils2=[15];
ndpt=[140];
AL1=generateL1mat(that1,that2,Apx,Apy,Apz,rangec,ncoils2(cc),16);

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
E2{des}=zeros(size(Emat2{1}));
for l=1:lmax
for m=1:1:l
E{des}=E{des}+coeff(globind(m,l)+stack)*Emat{globind(m,l)+stack};
E2{des}=E2{des}+coeff(globind(m,l)+stack)*Emat2{globind(m,l)+stack};
end
for m=0:1:l
E{des}=E{des}+coeff(globind(m,l))*Emat{globind(m,l)};
E2{des}=E2{des}+coeff(globind(m,l))*Emat2{globind(m,l)};
end
end


Im(des)=50/E2{des}(1,1,Ndep+1,1)
E{des}=sqrt(sum(abs(E{des}).^2,4));
E{des}=abs(Im(des))*E{des};
V(des)=sum(E{des}(:)>=50);
con(des)=max(E{des}(:));
Jv=Im(des)*Jcalc(:);
W(des)=Jv'*Lmat*Jv/2;
hold on
end
if des>=1
Wi{cc}=W;
Si{cc}=10^-2*V/(Ndep/10);
Si{cc}=Si{cc}(Wi{cc}~=0);
Wi{cc}=Wi{cc}(Wi{cc}~=0);
[Wi{cc} order]=sort(Wi{cc});
Si{cc}=Si{cc}(order)
plot(Wi{cc},Si{cc})

pause(5)
hold on
end

end
%%
close all
plot(Wi{cc},Si{cc},'o-')
