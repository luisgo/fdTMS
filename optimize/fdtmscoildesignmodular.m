
options=defaultoptimopts();

function [x, fval, exitflag, output]=runoptimisation(Aloo,jacloo,slackloo, ...
    Avector,Lmat,ncoils,stv,env,options,pererr,Nang,condi,Wmax,Jmax,Srmax,vmat,harmfile)

for loo=1:3
%% adaptive refinement    
    NA=nnz(slack==1);
    x=zeros([2*ncoils+NA(3) 1]);
    rangec=[];

%% cost function
jac=jac*10000;%volume at most 1
NA=nnz(slack==1);
f=zeros([ncoils+NA(3)+ncoils 1]);
f(ncoils+1:end-ncoils) = jac(slack==1);

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
intcon = ncoils+1:ncoils+NA;
%% Write the linear inequality constraints. 




%% write peak Jconstraints
%load  indevals.mat inde Jparty ASR;

load(harmfile,'Apx','Apy','Apz','that1','that2'); 
if numel(rangec)~=0
Aine1=generateL1mat(that1(thatrange,:),that2(thatrange,:),Apx,Apy,Apz,rangec,ncoils,16);
bine1(1:numel(Aine1(:,1)))=Jmax;
Aine1(1:neq*16,1:ncoils)=Aine1*vmat;
else
Aine1=[];
bine1=[];
end

%% Write the linear equality constraints. 
Aeq = zeros(1,ncoils+NA(3)+ncoils);
Aeq(1,1:ncoils) = -Avector(1,1:ncoils,1)*condi;
Aeq=sparse(Aeq);
beq(1)=-50.000;
Aeq(2,stvar+1:end)=1;
beq(2)=Wmax;
%Aeq(3,1:ncoils) = -Avector(1,1:ncoils,2)*condi;
%beq(3)=-99.99;
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

[x, fval, exitflag, output]  = intlinprog(f,intcon,sparse([Aine;Aine2]),sparse([bine;bine2]),Aeq,beq,lb,ub,xxx,options);
%%
if exitflag<=0
x=zeros([2*ncoils+NA(3) 1]);
break;
end
%%
end

initj=x(1:ncoils);
truequadratic=initj'*LL*initj
test=vmat*initj;
load(harmfile,'rangec')
Jmag=sqrt((Apx(:,1:ncoils)*test).^2+...
     (Apy(:,1:ncoils)*test).^2+...
     (Apz(:,1:ncoils)*test).^2);
 thatrange=Jmag(rangec)>Jmax;
 rangec=rangec(Jmag(rangec)>Jmax);
zslack=sum(x(stvar+1:end))
pause(2)

end











