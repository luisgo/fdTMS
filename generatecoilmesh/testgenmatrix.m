function testgenmatrix
x=.9;
N=5;
[P, dPx] = gen_matrix(N, x);
P2=zeros(N+1,N);
dPx2=zeros(N+1,N);

for i=0:N
    legendrep(i,N,x)
P2(i+1:(N+1),i+1)=legendrep(i,N,x);
dPx2(i+1:(N+1),i+1)=legendrepprime(i,N,x,legendrep(i,N,x))*(1-x.^2).^(-1/2);
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

function P=legendrep2(m,l,x)
%computes all legendre polynomials in an array P^{m}_{m,m+1,...,l} 0=<m<=l
P=zeros([(l-m)+1 length(x(:))]);
x=reshape(x(:),[1 length(x(:))]);
fac2=prod(1:2:(2*m-1));
P(1,:)=(-1)^(m)*fac2*(1-x.^2).^(m/2-1/2);
if m<l
P(2,:)=x*(2*m+1).*P(1,:);
end
for lval=m+2:l
P(lval-m+1,:)=((2*lval-1)*x.*P(lval-m,:)-(lval+m-1)*P(lval-m-1,:))/(lval-m);
end
end

function Pprime=legendrepprime(m,l,x,P)
%computes all legendre polynomials in an array P'^{m}_{m,m+1,...,l} 0=<m<=l
Pprime=zeros([(l-m)+1 length(x(:))]);
x=reshape(x(:),[1 length(x(:))]);
fac2=prod(1:2:(2*m-1));
if m~=0
Pprime(1,:)=m/2*(-1)^(m)*fac2*(-2*x).*(1-x.^2).^(m/2-1/2); 
else
Pprime(1,:)=0;     
end
if m<l
Pprime(2,:)=(2*m+1).*P(1,:).*(1-x.^2).^(1/2)+(2*m+1).*x.*Pprime(1,:);
end
for lval=m+2:l
Pprime(lval-m+1,:)=((2*lval-1)*P(lval-m,:).*(1-x.^2).^(1/2)+...
    (2*lval-1)*x.*Pprime(lval-m,:)-(lval+m-1)*Pprime(lval-m-1,:))/(lval-m);
end
Pprime=-Pprime;%-(1-x.^2)^1/2 is negative sin

end