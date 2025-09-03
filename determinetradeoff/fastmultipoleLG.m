function green=fastmultipoleLG(rs,rt,lmax)
%define source multipoles
rads=sqrt(sum(rs.^2,2));
phis= atan2(rs(:,2),rs(:,1));
thetas = acos(rs(:,3)./rads);
phis(rads==0)=0;
thetas(rads==0)=0;
[Ys]=Yvalconj(lmax,rads',thetas',phis');

%define target multipoles
radt=sqrt(sum(rt.^2,2));
phit= atan2(rt(:,2),rt(:,1));
thetat = acos(rt(:,3)./radt);
phit(radt==0)=0;
thetat(radt==0)=0;
[Yt]=Yval(lmax,radt',thetat',phit');
for i=1:3
mi(i)=min(cat(1,rs(:,i),rt(:,i)));
ma(i)=max(cat(1,rs(:,i),rt(:,i)));
end
del=(ma-mi);
green=zeros([numel(radt),1]);
green=fmmfunc(rs,Ys,rt,Yt,mi,ma,green);
for i=1:numel(green)
    green(i)=sum(Ys(:).*Yt(:,i));
end
end
function green=fmmfunc(rs,Ys,rt,Yt,mi,ma,del,green);

(rs(:,1)>mi(:))
mi+del

end
function [Y]=Yval(lmax,rad,theta,phi)
globind=@(l,m) (m+l)+1+l^2;
Y=zeros([(lmax+1)^2,length(theta(:))]);
rho=ones([lmax+1,numel(rad)]);
rho(1,:)=1./rad;
for l=1:lmax
rho(l+1,:)=rho(l,:)./rad;
end
for m=0:lmax
    legen=legendrep(m,lmax,cos(theta(:)));
    for innerloop=m:lmax
        l=innerloop;
        Const=sqrt(factorial(l-abs(m))/factorial(l+abs(m)));
        Y(globind(l,m),:)=Const*rho(l+1,:).*legen(l-abs(m)+1,:).*exp(1i*m*phi);
        Y(globind(l,-m),:)=Const*rho(l+1,:).*legen(l-abs(m)+1,:).*exp(-1i*m*phi);
    end
end

end


function [Y]=Yvalconj(lmax,rad,theta,phi)
globind=@(l,m) (m+l)+1+l^2;
Y=zeros([(lmax+1)^2,length(theta(:))]);

rho=ones([lmax+1,numel(rad)]);
for l=1:lmax
rho(l+1,:)=rho(l,:).*rad;
end
for m=0:lmax
    legen=legendrep(m,lmax,cos(theta(:)));
    for innerloop=m:lmax
        l=innerloop;
        Const=sqrt(factorial(l-abs(m))/factorial(l+abs(m)));
        Y(globind(l,-m),:)=Const*rho(l+1,:).*legen(l-abs(m)+1,:).*exp(1i*m*phi);
        Y(globind(l,m),:)=Const*rho(l+1,:).*legen(l-abs(m)+1,:).*exp(-1i*m*phi);
    end
end

end

function P=legendrep(m,l,x)
%computes all legendre polynomials in an array P^{m}_{m,m+1,...,l} 0=<m<=l
P=zeros([(l-m)+1 length(x(:))]);
x=reshape(x(:),[1 length(x(:))]);
fac2=prod(1:2:(2*m-1));
P(1,:)=(-1)^(m)*fac2*(1-x.^2).^(m/2); %P^m_l l=m
if m<l
P(2,:)=x.*(2*m+1).*P(1,:);%P^m_l l=m+1
end
for lval=m+2:l
P(lval-m+1,:)=((2*lval-1)*x.*P(lval-m,:)-(lval+m-1)*P(lval-m-1,:))/(lval-m); %P^m_l l=m+1+loopit
end

end

