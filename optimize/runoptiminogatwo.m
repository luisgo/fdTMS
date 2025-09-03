function ncoils = runoptiminoga(depth,energy,Jmax,Srmax,ncoils,harmfile,filename)
lmax=20;
filenamesav=strcat(filename,'d',num2str(depth(1)),'ncoils',num2str(ncoils));
spmin=codeupdaterefine(ncoils,depth,energy,Jmax,Srmax,harmfile,filenamesav);
end


function spmin=codeupdaterefine(ncoils,ndepth,Menergy,Jmax,Srmax,harmfile,filenamesav)
spmin=1;
%% inputs
stv=[40 45 45];
env=[70 65 65];
lmax=30;
condi=10^3;%parameter to keep curr value low
Nang=16;%number of directions
pererr=50*ncoils;%percent error in energy estimation
options=0;
%% setup
%[Aloo,Avector,jacloo,~,Lmat,slackloo]=planeintegratorrefine(ndepth,ncoils);
[Aloo,Avector,jacloo,~,Lmat,slackloo,vmat]=planeintegratorrefine2layer(ndepth,ncoils,harmfile);

startv=0;
%%
for jj=1:length(Menergy)
%% run optimisation
[x, fval, exitflag, output]=runoptimisation(Aloo,jacloo,slackloo, ...
    Avector,Lmat,ncoils,stv,env,options,pererr,Nang,condi,Menergy(jj),Jmax(jj)/condi,Srmax(jj)/condi,vmat,harmfile);
exflag(jj)=exitflag;
%% post-process
x=condi*x;

x(1:ncoils)=vmat*x(1:ncoils);
Jcalc=x;
Jcalc2{jj}=x(1:ncoils);
x2{jj}=x;
 save(strcat(filenamesav,'lmax',num2str(lmax),...
'.mat'),'-v7.3','Jcalc','Jcalc2','Menergy','Lmat','exflag');
end
end

function [x, fval, exitflag, output]=runoptimisation(Aloo,jacloo,slackloo, ...
    Avector,Lmat,ncoils,stv,env,options,pererr,Nang,condi,Wmax,Jmax,Srmax,vmat,harmfile)

rangec=[];

for loo=1:3
%% adaptive refinement    
if loo==1
    A=Aloo{loo};
    jac=jacloo{loo};
    slack=slackloo{loo};
    refl=ones([numel(jac) 1]);
    tind=1:numel(jac);
    NA(3)=nnz(slack==1);
    x=zeros([2*ncoils+NA(3) 1]);
else
    if loo==2
   [A,jac,slack,refl,tind]...
    =generatenewmesh(Aloo,jacloo,slackloo,refl,tind,A,x,ncoils,[stv(1) env(2)]/condi,[48 90]/condi);
    elseif loo==3
   [A,jac,slack,refl,tind]...
    =generatenewmesh(Aloo,jacloo,slackloo,refl,tind,A,x,ncoils,[stv(2) env(2)]/condi,[48 85]/condi);
    elseif loo==4
   [A,jac,slack,refl,tind]...
    =generatenewmesh(Aloo,jacloo,slackloo,refl,tind,A,x,ncoils,[stv(3) env(3)]/condi,[48 70]/condi);
    end
end


%% cost function
jac=jac*10000;%volume at most 1
NA(3)=nnz(slack==1);
f=zeros([ncoils+NA(3)+ncoils 1]);
f(ncoils+1:end-ncoils) = jac(slack==1);

%% write peak Jconstraints

load(harmfile,'Apx','Apy','Apz','that1','that2'); 
load  indevals.mat inde Jparty ASR;
clear bine2;
if numel(rangec)~=0
    
    AL1=generateL1mat(that1(thatran,:),that2(thatran,:),Apx,Apy,Apz,rangec,ncoils,16);
    neq=numel(rangec);
Aine2=zeros([neq*16 ncoils+NA(3)+ncoils]);
%Aine2(1:numel(inde),1:ncoils)=Jparty(:,1:ncoils)*vmat;
%bine2(1:numel(inde))=Srmax;
bine2(1:neq*16)=Jmax;
Aine2(1:neq*16,1:ncoils)=AL1*vmat;
else
Aine2=[];
bine2=[];
end
%% write bounds

stvar=ncoils+NA(3);
lb = zeros(ncoils+NA(3)+ncoils,1);
ub = ones(ncoils+NA(3)+ncoils,1);
for i=1:ncoils
lb(i)=-10^5;
ub(i)=10^5;
end
LL=Lmat*condi^2/2;
for i=1:ncoils
lb(i)=-sqrt(Wmax/LL(i,i));
ub(i)=sqrt(Wmax/LL(i,i));
end
lb(stvar+1:end)=0;
ub(stvar+1:end)=Wmax;
%% integer constraints
intcon = ncoils+1:ncoils+NA(3);
%% Write the linear inequality constraints. 
thi=0:2*pi/Nang:2*pi-2*pi/Nang;
nsamp=round(sqrt(4*pererr));
Aine=zeros([Nang*numel(slack)+nsamp*ncoils  ncoils+NA(3)+ncoils]);
%write field constraints
ct=1;
bine=50*ones([Nang*numel(slack)+nsamp*ncoils 1]);
slackid=zeros([Nang*NA(3) 1]);
for i=1:numel(slack)
    th=thi+rand(1)*2*pi/Nang;%shift it so points dont line up
for j=1:Nang
    Aine(j+(i-1)*Nang,1:ncoils)=(A(1,:,i)*sin(th(j))+A(2,:,i)*cos(th(j)))*condi;
        if slack(i)==1
    Aine(j+(i-1)*Nang,ncoils+ct)=-50;
    slackid(j+(ct-1)*Nang)=j+(i-1)*Nang;
        end
if slack(i)==2
bine(j+(i-1)*Nang)=100;
end
end
if slack(i)==1
ct=ct+1;
end
end
Aine=sparse(Aine);
%energy constraints
steqn=Nang*numel(slack);
for i=1:ncoils
    Xmax=ub(i);
    Xmin=lb(i);
    for j=1:nsamp
        lam=(j-1/2)/nsamp;
        xoi=Xmin*lam+Xmax*(1-lam);
    Aine(steqn+j+(i-1)*nsamp,i)=2*LL(i,i)*xoi;
    Aine(steqn+j+(i-1)*nsamp,stvar+i)=-1;
    bine(steqn+j+(i-1)*nsamp)=LL(i,i)*xoi^2;
    end
end

%% Write the linear equality constraints. 
Aeq = zeros(1,ncoils+NA(3)+ncoils);
Aeq(1,1:ncoils) = -Avector(1,1:ncoils,1)*condi;
Aeq=sparse(Aeq);
beq(1)=-50.000;
Aeq(2,stvar+1:end)=1;
beq(2)=Wmax;
Aeq(3,1:ncoils) = -Avector(1,1:ncoils,2)*condi;
beq(3)=-100/sqrt(2);
%% run optimization
%%

yo=bine-Aine(:,1:ncoils)*x(1:ncoils);
if nnz(x)==0
xxx=[];
else
xxx=x(1:ncoils);
for i=1:NA(3)
if nnz(yo(slackid((i-1)*Nang+1:i*Nang))<0)~=0
xxx(ncoils+i)=1;
end
end
xxx(ncoils+NA(3)+1:2*ncoils+NA(3))=x(end-ncoils+1:end);
conscheck=nnz((bine(:)-Aine*xxx(:))<-10^-10)
conscheckeq=sum((beq(:)-Aeq*xxx(:)))
initcost=f'*xxx(:)

end
if loo==4
save res.mat
end
options=optimoptions('intlinprog');
options.AbsoluteGapTolerance= 0;
options.BranchRule= 'maxpscost';
options.ConstraintTolerance= 1.0000e-04;
options.CutGeneration= 'advanced';
 options.CutMaxIterations= 25;
options.Display= 'iter';
options.Heuristics= 'rins-diving';
options.HeuristicsMaxNodes= 500;
options.IntegerPreprocess= 'advanced';
options.IntegerTolerance= 1.0000e-05;
options.LPMaxIterations= 10^6;
options.LPOptimalityTolerance= 1.0000e-07;
options.MaxFeasiblePoints= Inf;
 options.MaxNodes= 10000000000;
 options.MaxTime= 6*7200;
options.NodeSelection= 'simplebestproj';
options.ObjectiveCutOff= Inf;
options.ObjectiveImprovementThreshold= 1.0000e-04;
options.OutputFcn= [];
options.PlotFcn= [];
options.RelativeGapTolerance= 1.0000e-04;
options.RootLPAlgorithm= 'dual-simplex';
 options.RootLPMaxIterations=10^6; 

[x, fval, exitflag, output]  = intlinprog(f,intcon,sparse([Aine;Aine2]),sparse([bine(:);bine2(:)]),Aeq,beq,lb,ub,xxx,options);
%%
if exitflag<=0
x=zeros([2*ncoils+NA(3) 1]);
break;
end
%%

initj=x(1:ncoils);
test=vmat*initj;
load(harmfile,'rangec')
Jmag=sqrt((Apx(:,1:ncoils)*test).^2+...
     (Apy(:,1:ncoils)*test).^2+...
     (Apz(:,1:ncoils)*test).^2);
 if loo~=1
     err=nnz(Jmag(rangec)>Jmax)
 thatran=(Jmag(rangec)>Jmax+thatran)>0;
 else
     err0=nnz(Jmag(rangec)>Jmax)
    thatran=(Jmag(rangec)>Jmax);  
 end
 rangec=rangec(thatran);
 
end

initj=x(1:ncoils);
truequadratic=initj'*LL*initj
test=vmat*initj;
load(harmfile,'rangec')
Jmag=sqrt((Apx(:,1:ncoils)*test).^2+...
     (Apy(:,1:ncoils)*test).^2+...
     (Apz(:,1:ncoils)*test).^2);
 if loo~=1
     err=nnz(Jmag(rangec)>Jmax)
 thatran=(Jmag(rangec)>Jmax+thatran)>0;
 else
     err0=nnz(Jmag(rangec)>Jmax)
    thatran=(Jmag(rangec)>Jmax);  
 end
 rangec=rangec(thatran);
 
zslack=sum(x(stvar+1:end))
pause(2)

end

function [A,Avector,jac,harmcoeff,Lmat,slack,vmat]=planeintegratorrefine2layer(ndepth,nkeep,harmfile)
load(harmfile,'Lmat2','harmcoeff','bigR','lmax');

ncoils=nkeep;
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
%get inductance
Lmat2=Lmat2(1:nkeep,1:nkeep);
[~,~,vmat]=svd(Lmat2);
Lmat=vmat'*Lmat2*vmat;
%%%choose sample points
loomax=4;
for loo=1:loomax
if loo==1
load baricooroNEW1.mat
elseif loo==2
load baricooroNEW2.mat
elseif loo==3
load baricooroNEW3.mat
elseif loo==4
load baricooroNEW4.mat
end
%get fields
matsize(loo)=numel(Zeval);
Bfield{loo}=Eevaluator(Xeval,Yeval,Zeval,bigR,lmax);%volume field
A{loo}=zeros([3 nkeep numel(Xeval)]); %volume mat
jactemp{loo}=jac;
slacktemp{loo}=slack;
end
jac=jactemp;
slack=slacktemp;
ct=1;
Zeval(1)=.07-0.0001*ndepth(1);
Zeval(2)=.07-0.0001*ndepth(2);
Xeval=zeros(size(Zeval));
Yeval=zeros(size(Zeval));
matsize2=numel(Zeval);
[Bfield2]=Eevaluator(Xeval,Yeval,Zeval,bigR,lmax); %surface field
Avector=zeros([3 ncoils numel(Xeval)]); %surface mat

stack=globind(lmax,lmax);
for icoils=1:ncoils
    coeff=zeros(size(harmcoeff{1}));
    for vloo=1:numel(vmat(:,1))
coeff=coeff+real(harmcoeff{vloo})*vmat(vloo,icoils);
    end
for l=1:lmax%populate field matrices
    for m=1:1:l
        for loo=1:loomax
        A{loo}(1:3,icoils,:)=A{loo}(1:3,icoils,:)...
            +permute(reshape(Bfield{loo}{globind(m,l)+stack},[matsize(loo) 1 3])*coeff(globind(m,l)+stack),[3 2 1]);
        end
        Avector(1:3,icoils,:)=Avector(1:3,icoils,:)...
            +permute(reshape(Bfield2{globind(m,l)+stack},[matsize2 1 3])*coeff(globind(m,l)+stack),[3 2 1]); 
    end
    for m=0:1:l
        for loo=1:loomax
        A{loo}(1:3,icoils,:)=A{loo}(1:3,icoils,:)...
            +permute(reshape(Bfield{loo}{globind(m,l)},[matsize(loo) 1 3])*coeff(globind(m,l)),[3 2 1]);
        end
        Avector(1:3,icoils,:)=Avector(1:3,icoils,:)...
            +permute(reshape(Bfield2{globind(m,l)},[matsize2 1 3])*coeff(globind(m,l)),[3 2 1]);
    end
end

end

end





function [Efield,Srpos,Srneg]=Eevaluator(Xeval,Yeval,Zeval,bigR,lmax)
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
stack=globind(lmax,lmax);
%change to spherical
nx=size(Xeval);
if numel(nx)<3;
nx(end+1:3)=1;
end
Reval=sqrt(Xeval(:).^2+Yeval(:).^2+Zeval(:).^2);
rad=Reval;
phi = atan2(Yeval(:),Xeval(:));
theta = acos(Zeval(:)./sqrt(Xeval(:).^2+Yeval(:).^2+Zeval(:).^2));

[Ypos,Yneg,Srpos,Srneg]=Yfunc2(lmax,theta,phi,rad);
Npt=numel(theta);

for l=0:lmax
    cons=-1/(2*l+1)/bigR^l;
    radrr(1,:,1)=cons*reshape(rad(:).^(l),[1,numel(rad),1]);
    Br=radrr.*Ypos{2}(globind(0,l),:,1);Bth=radrr.*Ypos{2}(globind(0,l),:,2);Bph=radrr.*Ypos{2}(globind(0,l),:,3);
        Efield{globind(0,l)}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(0,l)}(:,:,:,1)=reshape(Bth(:),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)}(:,:,:,2)=reshape(Bph(:),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)}(:,:,:,3)=0;
        
        Efield{globind(0,l)+stack}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(0,l)+stack}(:,:,:,1)=reshape(Bth(:),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)+stack}(:,:,:,2)=reshape(Bph(:),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(0,l)+stack}(:,:,:,3)=0;
    for m=1:l
    Bth=radrr.*Ypos{2}(globind(m,l),:,2);Bph=radrr.*Ypos{2}(globind(m,l),:,3);
        Efield{globind(m,l)}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(m,l)}(:,:,:,1)=reshape(Bth(:),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)}(:,:,:,2)=reshape(Bph(:),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)}(:,:,:,3)=0;    
        
    Bth=radrr.*Yneg{2}(globind(m,l),:,2);Bph=radrr.*Yneg{2}(globind(m,l),:,3);
        Efield{globind(m,l)+stack}=zeros([nx(1) nx(2) nx(3) 3]);
        Efield{globind(m,l)+stack}(:,:,:,1)=reshape(Bth(:),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)+stack}(:,:,:,2)=reshape(Bph(:),[nx(1) nx(2) nx(3) 1]);
        Efield{globind(m,l)+stack}(:,:,:,3)=0;  
    end
end
end
function [Ypos,Yneg,Srpos,Srneg]=Yfunc2(lmax,theta,phi,rad)
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
Npt=length(theta(:));
theta=reshape(theta(:),[1 Npt]);
phi=reshape(phi(:),[1 Npt]);
rad=reshape(rad(:),[1 Npt]);
Ypos{2}=zeros([globind(lmax,lmax),Npt 3]);
Yneg{2}=zeros([globind(lmax,lmax),Npt 3]);
ct=1;
for m=1:lmax
Pvals=legendrep(m,lmax,cos(theta));
Pvals2=legendrep2(m,lmax,cos(theta)); %ensure no 0/0 at a small cost
Pprime=legendrepprime(m,lmax,cos(theta),Pvals);%this prime has the -sin(theta) integrated
for inloop=m:lmax
    l=inloop;
    harmind=globind(m,l);
    locind=inloop-m+1;
    Const=1/sqrt(l*(l+1))*sqrt(((2*l+1)/(2*pi)*factorial(l-m)/factorial(l+m)));    
Ypos{2}(harmind,:,2)=-Const*Pvals2(locind,:).*(-m*sin(m*phi)); %derivative w.r.t. phi
Ypos{2}(harmind,:,3)=(Const).*Pprime(locind,:).*cos(m*phi); %derivative w.r.t. theta
Yneg{2}(harmind,:,2)=Const*Pvals2(locind,:).*(m*cos(m*phi)); %derivative w.r.t. phi
Yneg{2}(harmind,:,3)=-(Const).*Pprime(locind,:).*sin(m*phi); %derivative w.r.t. theta

end
end
Pvals=legendrep(0,lmax,cos(theta));
Pprime=legendrepprime(0,lmax,cos(theta),Pvals);
for inloop=2:lmax+1
    l=inloop-1;
    harmind=globind(0,l);m=0;
    Const=1/sqrt(l*(l+1))*sqrt(((2*l+1)/(4*pi)*factorial(l-m)/factorial(l+m)));
Ypos{2}(harmind,:,3)=Const*Pprime(inloop,:);
Yneg{2}(harmind,:,3)=Ypos{2}(harmind,:,3);
end
[Srpos,Srneg]=Yval(lmax,theta,phi);
end
function [Ypos,Yneg]=Yval(lmax,theta,phi)
globind=@(m,l) l+1+m*(lmax-(m-1)/2);
Ypos=zeros([globind(lmax,lmax),length(theta(:))]);
Yneg=zeros([globind(lmax,lmax),length(theta(:))]);
for m=1:lmax
    legen=legendrep(m,lmax,cos(theta(:)));
    for innerloop=m:lmax
        l=innerloop;
        harmind=globind(m,l);
        Const=-1/sqrt(l*(l+1))*sqrt(((2*l+1)/(2*pi)*factorial(l-m)/factorial(l+m)));
        Ypos(harmind,:)=Const*legen(l-m+1,:).*cos(m*phi);
        Yneg(harmind,:)=-Const*legen(l-m+1,:).*sin(m*phi);
    end
end
legen=legendrep(0,lmax,cos(theta(:)));
    for innerloop=1:lmax
        l=innerloop;
        m=0;
        Const=-1/sqrt(l*(l+1))*sqrt(((2*l+1)/(2*pi)*factorial(l-m)/factorial(l+m)));
Ypos(l+1,:)=Const*legen(l+1,:)/sqrt(2);
Yneg(l+1,:)=Ypos(innerloop+1,:);
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