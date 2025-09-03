function [Wi,Si,peakJ,peakSR,Sd,Wd]=anycoildeterminetradeoff_f(harmfile,designfile,Ndep,ncoils);
load(harmfile);
if  exist('shortcut.mat')~=0
load shortcut.mat Emat2 Emat bigR2;
else
[Yeval,Xeval,Zeval]=meshgrid(-.07:.001:.07,-.07:.001:.07,0.0305:.001:.07);
[Emat]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
[Yeval,Xeval,Zeval]=meshgrid(0,0,.07:-.0001:.01);
[Emat2]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
bigR2=bigR;
save shortcut.mat -v7.3 Emat2 Emat bigR2;

end
if bigR2~=bigR 
[Yeval,Xeval,Zeval]=meshgrid(-.07:.001:.07,-.07:.001:.07,0.0305:.001:.07);
[Emat]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
[Yeval,Xeval,Zeval]=meshgrid(0,0,.07:-.0001:.01);
[Emat2]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
bigR2=bigR;
save shortcut.mat -v7.3 Emat2 Emat bigR2;
end
[Yeval,Xeval,Zeval]=meshgrid(-.07:.001:.07,-.07:.001:.07,0.0305:.001:.07);
reval=sqrt(Xeval.^2+Yeval.^2+Zeval.^2)<=.07;
clear Xeval Yeval Zeval;
for i=1:numel(Emat)
Emat{i}(:,:,:,1)=Emat{i}(:,:,:,1).*reval;
Emat{i}(:,:,:,2)=Emat{i}(:,:,:,2).*reval;
Emat{i}(:,:,:,3)=Emat{i}(:,:,:,3).*reval;
end
AL1=generateL1mat(that1,that2,Apx,Apy,Apz,rangec,ncoils,16);

    clear W V
load(designfile)
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


Im(des)=50/E2{des}(1,1,Ndep+1,1);
E2{des}=Im(des)*E2{des};

E{des}=sqrt(sum(abs(E{des}).^2,4));
E{des}=abs(Im(des))*E{des};
V(des)=sum(E{des}(:)>=50);
con(des)=max(E{des}(:));
con2(des)=max(E2{des}(1,1,1,1));
con
con2
Jv=Im(des)*Jcalc(:);
peakSR=max(abs(Asym*Atrans(:,1:ncoils)*Jcalc(:)));
peakJ(des)=max(AL1*Jcalc(:));
W(des)=Jv'*Lmat*Jv/2;
Sd{des}=zeros([Ndep 1]);
Wd{des}=zeros([Ndep 1]);
for ides=1:Ndep
Sd{des}(ides)=10^-1*sum(E{des}(:)>=E2{des}(1,1,ides+1,1))/ides;
Wd{des}(ides)=W(des)*(50/E2{des}(1,1,ides+1,1))^2;
end
hold on
end
if des>=1
Wi=W;
Si=10^-2*V/(Ndep/10);
Si=Si(Wi~=0);
Wi=Wi(Wi~=0);
[Wi order]=sort(Wi);
Si=Si(order)
plot(Wi,Si)

pause(5)
hold on
end

end