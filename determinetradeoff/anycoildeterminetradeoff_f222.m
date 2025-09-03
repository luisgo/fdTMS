function [Wi,Si,peakJ,peakSR,peakEscalp,EscalpS12,...
    Ndep,Evol,Edecay,Escalp,Sisq,Ndepsq,Sd,Wd]=anycoildeterminetradeoff_f(harmfile,designfile,Ndep,ncoils)
load(harmfile);

load shortcut.mat jac3 Emat3 Emat2 Emat bigR2;

 
if bigR2~=bigR
[Yeval,Xeval,Zeval]=meshgrid(-.07:.001:.07,-.07:.001:.07,0.0305:.001:.07);
[Emat]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
[Yeval,Xeval,Zeval]=meshgrid(0,0,.07:-.0001:.01);
[Emat2]=Eevaluator2(Xeval,Yeval,Zeval,bigR,lmax);
bigR2=bigR;

load('C:\Users\ljg24\Desktop\Coiloptcode\optimize\baricooroNEW4.mat');
jac3=jac/sum(jac)*(4*pi*8.5^2);
[Emat3]=Eevaluator2(Xeval*85/70,Yeval*85/70,Zeval*85/70,bigR,lmax);
save shortcut.mat -v7.3 jac3 Emat3 Emat2 Emat bigR2;
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
E3{des}=zeros(size(Emat3{1}));
for l=1:lmax
for m=1:1:l
E{des}=E{des}+coeff(globind(m,l)+stack)*Emat{globind(m,l)+stack};
E2{des}=E2{des}+coeff(globind(m,l)+stack)*Emat2{globind(m,l)+stack};
E3{des}=E3{des}+coeff(globind(m,l)+stack)*Emat3{globind(m,l)+stack};
end
for m=0:1:l
E{des}=E{des}+coeff(globind(m,l))*Emat{globind(m,l)};
E2{des}=E2{des}+coeff(globind(m,l))*Emat2{globind(m,l)};
E3{des}=E3{des}+coeff(globind(m,l))*Emat3{globind(m,l)};
end
end

% if Ndep==96
%     Ndep=47;
% elseif Ndep==125
%     Ndep=62;
% elseif Ndep==128
%     Ndep=63;
% elseif Ndep==142
%     Ndep=71;
% elseif Ndep==168
%     Ndep=86;
% elseif Ndep==199
%     Ndep=102;
% end
% 
% if Ndep==47
%     Ndep=96;
% elseif Ndep==62
%     Ndep=125;
% elseif Ndep==63
%     Ndep=128;
% elseif Ndep==71
%     Ndep=142;
% elseif Ndep==86
%     Ndep=168;
% elseif Ndep==102
%     Ndep=199;
% end

Im(des)=100/2/E2{des}(1,1,Ndep+1,1)
E{des}=sqrt(sum(abs(E{des}).^2,4));
E{des}=abs(Im(des))*E{des};
E2{des}=abs(Im(des))*E2{des};
E3{des}=abs(Im(des))*E3{des};
V(des)=sum(E{des}(:)>=100/2);
Vsq(des)=sum(E{des}(:)>=100/sqrt(2));
con(des)=max(E{des}(:));
con2(des)=max(E2{des}(1,1,1,1));
Evol{des}=E{des};
Edecay{des}=E2{des};
Escalp{des}=E3{des};
Ndepsq(des)=nnz(Edecay{des}(1,1,:,1)>100/sqrt(2))-1;
con
con2
Jv=Im(des)*Jcalc(:);
peakJ(des)=max(AL1*Jcalc(:));
peakSR(des)=max(abs(Asym*Atrans(:,1:ncoils)*Jcalc(:)));

E3{des}=Im(des)*sqrt(sum(abs(E3{des}).^2,4));
load barisamp4.mat
trisurf(t2p,p(:,1),p(:,2),p(:,3),E3{des})
peakEscalp(des)=max(E3{des});
EscalpS12(des)=sum(jac3(:).*(E3{des}(:)>=peakEscalp(des)/2))*4;
W(des)=Jv'*Lmat*Jv/2;
hold on

Sd{des}=zeros([Ndep 1]);
Wd{des}=zeros([Ndep 1]);
for ides=1:Ndep
Sd{des}(ides)=10^-1*sum(E{des}(:)>=E2{des}(1,1,ides+1,1))/ides;
Wd{des}(ides)=W(des)*(50/E2{des}(1,1,ides+1,1))^2;
end

end
if des>=1
Wi=W;
Si=10^-2*V/(Ndep/10);
Si=Si(Wi~=0);
Sisq=10^-2*Vsq./(Ndepsq/10);
Sisq=Sisq(Wi~=0);
Wi=Wi(Wi~=0);
[Wi order]=sort(Wi);
Si=Si(order)
Sisq=Sisq(order);
Ndepsq=Ndepsq(order);
plot(Wi,Si)

pause(5)
hold on
end

end