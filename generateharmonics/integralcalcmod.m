function [A]=integralcalcmod(tvols,nhats,pop,psp,tol)
A=zeros([3 3]);

s(:,1)=[1;1;-1;-1];
s(:,2)=[1;-1;1;-1];
s(:,3)=[-1;-1;1;1];

po=pop;ps=psp;
for ko=1:4
        po(:,1)=s(ko,1)*pop(:,1);
        po(:,2)=s(ko,2)*pop(:,2);       
for ks=1:4
        ps(:,1)=s(ks,1)*psp(:,1);
        ps(:,2)=s(ks,2)*psp(:,2);

ceno=sum(po,1)/3;cens=sum(ps,1)/3;
rad=norm(ceno-cens);
if rad<3*tol
    nquad=4;
    a=3/5;b=1/5;
    qpt(1,:)=[a b b];
    qpt(2,:)=[b a b];
    qpt(3,:)=[b b a];
    qpt(4,:)=[1/3 1/3 1/3];
    qwt(1:3)=25/96*2;
    qwt(4)=-9/32*2;
else
    nquad=1;
    qpt=[1/3 1/3 1/3];
    qwt=1;    
end
if rad>tol
    for qo=1:nquad
        ceno=po(1,:)*qpt(qo,1)+po(2,:)*qpt(qo,2)+po(3,:)*qpt(qo,3);
    for qs=1:nquad
        cens=ps(1,:)*qpt(qs,1)+ps(2,:)*qpt(qs,2)+ps(3,:)*qpt(qs,3);
        rad=norm(ceno-cens);
    for is=1:3
        basisint=s(ks,3)*(cens-ps(is,:));
        for io=1:3 
            A(io,is)=A(io,is)+s(ko,3)*qwt(qo)*qwt(qs)*sum((ceno-po(io,:)).*basisint)/rad;
        end
    end
    end
    end
else
for qo=1:nquad
ceno=po(1,:)*qpt(qo,1)+po(2,:)*qpt(qo,2)+po(3,:)*qpt(qo,3);    
[K1,boldk2]=pasiintegration(ceno,ps,-1);
for is=1:3
basisint=2*boldk2{is}(2,:);
for io=1:3
    A(io,is)=A(io,is)+s(ko,3)*s(ks,3)*qwt(qo)*sum((ceno-po(io,:)).*basisint);
end
end
end
end
end
end
A=A/4;
end

function [K1,boldk2]=pasiintegration(ro,q,nmax)

ro=reshape(ro,[1 3]);
cen=sum(q,1)/3;
dir(3,:)=q(2,:)-q(1,:);
dir(2,:)=q(1,:)-q(3,:);
dir(1,:)=q(3,:)-q(2,:);

%l3 l2 l1
l(3)=norm(dir(3,:));l(2)=norm(dir(2,:));l(1)=norm(dir(1,:));
%%% local coordinate system calc (uhat,vhat,nhat)
nhat=reshape([dir(2,2)*dir(3,3)-dir(2,3)*dir(3,2) ...
              dir(2,3)*dir(3,1)-dir(2,1)*dir(3,3) ...
              dir(2,1)*dir(3,2)-dir(2,2)*dir(3,1)],[1 3]);
A2=norm(nhat);
nhat=nhat/A2;
uhat=dir(3,:)/l(3);
vhat=reshape([nhat(2)*uhat(3)-nhat(3)*uhat(2) ...
              nhat(3)*uhat(1)-nhat(1)*uhat(3) ...
              nhat(1)*uhat(2)-nhat(2)*uhat(1)],[1 3]);
  %%%
u3=-sum(dir(2,:).*dir(3,:))/l(3);
v3=A2/l(3);
%%%zeros calcs
u0=sum((ro-q(1,:)).*uhat);
v0=sum((ro-q(1,:)).*vhat);
w0=sum((ro-q(1,:)).*nhat);
%%%s calcs
s_m(1)=-((l(3)-u0)*(l(3)-u3)+v0*v3)/l(1);
s_m(2)=-(u3*(u3-u0)+v3*(v3-v0))/l(2);
s_m(3)=-u0;
s_p=s_m+l;
t0(1)=(v0*(u3-l(3))+v3*(l(3)-u0))/l(1);
t0(2)=(u0*v3-v0*u3)/l(2);
t0(3)=v0;
R0=sqrt(t0.^2+w0^2);
Rp=[norm(ro-q(3,:)) norm(ro-q(1,:)) norm(ro-q(2,:))];
Rm=[Rp(3) Rp(1) Rp(2)];
rho=ro-w0*nhat;
%%%%%
m=zeros(3);
for i=1:3
m(i,:)=(q(mod(i,3)+1,:)-cen);
m(i,:)=m(i,:)-sum(m(i,:).*dir(i,:))*dir(i,:)/l(i)^2;
m(i,:)=m(i,:)/norm(m(i,:));
end

K1=zeros([(nmax+1)/2+2 1]);
nplus2arr=(1:2:nmax+2);
narr=(-1:2:nmax+2);
for i=1:3
I{i}=Iord(Rp(i),s_p(i),Rm(i),s_m(i),R0(i),nmax+2);
end
if w0==0
    K1(1)=0;%;-sum(s_p./(t0.*Rp)-s_m./(t0.*Rm));
    K1(2:end)=1./nplus2arr(:).*(t0(1)*I{1}(1:end-1)+t0(2)*I{2}(1:end-1)+t0(3)*I{3}(1:end-1));
else
    for j=1:(nmax+1)/2+2
        if j==1;
            beta=atan((t0.*s_p)./(R0.^2+abs(w0)*Rp))-...
                 atan((t0.*s_m)./(R0.^2+abs(w0)*Rm));
             K1(1)=1/abs(w0)*sum(beta);
        else
             K1(j)=1/nplus2arr(j-1)*((narr(j-1))*w0^2*K1(j-1)+...
                 (t0(1)*I{1}(j-1)+t0(2)*I{2}(j-1)+t0(3)*I{3}(j-1)));
        end
    end
end
for i=1:3
boldI(:,i)=m(1,i)*I{1}+m(2,i)*I{2}+m(3,i)*I{3};
end
for i=1:3
    vec=-ro+q(i,:);
    vec2=(rho-q(i,:));
for j=1:(nmax+1)/2+2
    boldk2{i}(j,:)=boldI(j,:)/narr(j)+vec2*K1(j);
    
end
end
K1=2*K1/A2;%used for charge
for i=1:3
boldk2{i}=boldk2{i}/A2;
end

end


function I=Iord(Rp,sp,Rm,sm,R0,nmax,t0,mint)
I=zeros([(nmax+1)/2+1 1]);
ct=1;
for n=-1:2:nmax
if n==-1
    I(ct,1)=(log((Rp+sp)/(Rm+sm)));ct=ct+1;
     if isfinite(I(ct-1))==0
         I(ct-1)=0;
     end
else
I(ct,1)=1/(n+1)*(sp*Rp^n-sm*Rm^n+n*R0^2*I(ct-1));ct=ct+1;
end
end
end
