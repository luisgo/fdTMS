function Doublelayer2D()
sigma=10;
N=1000;
p(:,1)=cos(0:2*pi/N:2*pi-2*pi/N);
p(:,2)=sin(0:2*pi/N:2*pi-2*pi/N);
ed(:,1)=1:N;
ed(:,2)=ed(:,1)+1;
ed(end,2)=ed(1,1);
ed=ed(end:-1:1,:);
ned=numel(ed)/2;
np=numel(p)/2;
A=zeros([ned,ned]);
A2=zeros([ned,ned]);
ro=(p(ed(:,1),:)+p(ed(:,2),:))/2;
plot(p(:,1),p(:,2))
hold on
scatter(ro(:,1),ro(:,2))
pause
for jv=1:ned   
for iv=1:ned

    if iv~=jv
        A(iv,jv)=-gradgreenfunc(ro(iv,:),p(ed(jv,1),:),p(ed(jv,2),:))/(2*pi);
        A2(iv,jv)=-gradgreenfunc2(ro(iv,:),p(ed(jv,1),:),p(ed(jv,2),:))/(2*pi);
    else
       A(iv,jv)=1/2;
       A2(iv,jv)=1/2;
    end
end
end
b=zeros([ned 1]);
ps=[1/2,0];
for iv=1:ned
b(iv)=-log(norm(ro(iv,:)-ps(1,:)))/(2*pi);
end
x=A\b;
close all
plot(x)

end

    function I=gradgreenfunc(ro,pm,pp)
        lhat=pp-pm;
        lhat=lhat/norm(lhat);
        uhat=[lhat(2),-lhat(1)];
        P0=(ro(1)-pm(1)).*uhat(1)+(ro(2)-pm(2)).*uhat(2);
        lp=(ro(1)-pp(1)).*lhat(1)+(ro(2)-pp(2)).*lhat(2);
        lm=(ro(1)-pm(1)).*lhat(1)+(ro(2)-pm(2)).*lhat(2);
        I=-(atan(lp/P0)-atan(lm/P0));
    end


    function I=gradgreenfunc2(ro,pm,pp)
        lhat=pp-pm;
        len=norm(lhat);
        lhat=lhat/len;
        uhat=[lhat(2),-lhat(1)];
        P0=(ro(1)-pm(1)).*uhat(1)+(ro(2)-pm(2)).*uhat(2);
        I=0;
        nquad=4;
         [x,qwt]=lgwt(nquad,0,1);
         qpt(:,1)=x(:);qpt(:,2)=1-x(:);
         qwt=qwt*len;
        for i=1:nquad
        pmid=pp*qpt(i,1)+pm*qpt(i,2);
        lav=(ro(1)-pmid(1)).*lhat(1)+(ro(2)-pmid(2)).*lhat(2);
        I=I+qwt(i)*P0/(lav^2+P0^2);
        end
    end
    
    
    function [x,w]=lgwt(N,a,b)
N=N-1;
N1=N+1; N2=N+2;

xu=linspace(-1,1,N1)';
y=cos((2*(0:N)'+1)*pi/(2*N+2))+(0.27/N1)*sin(pi*xu*N/N2);
L=zeros(N1,N2);
Lp=zeros(N1,N2);
y0=2;
while max(abs(y-y0))>eps
    L(:,1)=1;
    Lp(:,1)=0;
    L(:,2)=y;
    Lp(:,2)=1;
    for k=2:N1
        L(:,k+1)=( (2*k-1)*y.*L(:,k)-(k-1)*L(:,k-1) )/k32;
    end
    Lp=(N2)*( L(:,N1)-y.*L(:,N2) )./(1-y.^2);   
    y0=y;
    y=y0-L(:,N2)./Lp;
end
x=(a*(1-y)+b*(1+y))/2;      
w=(b-a)./((1-y.^2).*Lp.^2)*(N2/N1)^2;
end