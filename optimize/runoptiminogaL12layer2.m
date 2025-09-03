function ncoils = runoptiminogaL12layer1(depth,energy,Jmax,gap,ncoils,harmfile,filename)
filenamesav=strcat(filename,'d',num2str(depth),'ncoils',num2str(ncoils));
spmin=codeupdaterefine(ncoils,depth,energy,Jmax,gap,harmfile,filenamesav);
end


function spmin=codeupdaterefine(ncoils,ndepth,Menergy,Jmax,gap,harmfile,filenamesav)
spmin=1;
%% inputs
stv=[40 45 45];
env=[70 65 65];
load(harmfile,'lmax');
condi=10^3;%parameter to keep curr value low
Nang=16;%number of directions
pererr=50*ncoils;%percent error in energy estimation
options=0;
%% setup
[Aloo,Avector,jacloo,~,Lmat,slackloo,vmat]=planeintegratorrefine2layer(ndepth,ncoils,harmfile);
startv=0;
%%
for jj=1:length(Menergy)
%% run optimisation
AL1=0;
[x, fval, exitflag, output]=runoptimisation(Aloo,jacloo,slackloo, ...
    Avector,Lmat,AL1,vmat,ncoils,stv,env,options,pererr,Nang,condi,Menergy(jj),Jmax(jj),gap(jj),harmfile);
exflag(jj)=exitflag;
%% post-process
x=condi*x;
x(1:ncoils)=vmat*x(1:ncoils);
Jcalc=x;
Jcalc2{jj}=x(1:ncoils);
x2{jj}=x;
 save(strcat(filenamesav,num2str(lmax),...
'ncoilstwolayer.mat'),'-v7.3','Jcalc','Jcalc2','Menergy','gap','Jmax','Lmat','exflag');
end
end

function [x, fval, exitflag, output]=runoptimisation(Aloo,jacloo,slackloo, ...
    Avector,Lmat,AL2,vmat,ncoils,stv,env,options,pererr,Nang,condi,Wmax,Jmax,gap,harmfile)
load  indevals.mat inde Jparty ASR;

Jmax=Jmax/condi;

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
load(harmfile,'that1','that2','Apx','Apy','Apz','rangec'); 

    boolvar=zeros([numel(that1(:,1)) 1]);
    ptt=1:numel(that1(:,1));
    ptt=ptt(:);
else
    if loo==2
   [A,jac,slack,refl,tind]...
    =generatenewmesh(Aloo,jacloo,slackloo,refl,tind,A,x(1:2*ncoils+NA(3)),ncoils,[stv(1) env(2)]/condi,[48 90]/condi);
    elseif loo==3
   [A,jac,slack,refl,tind]...
    =generatenewmesh(Aloo,jacloo,slackloo,refl,tind,A,x(1:2*ncoils+NA(3)),ncoils,[stv(2) env(2)]/condi,[48 85]/condi);
    elseif loo==4
   [A,jac,slack,refl,tind]...
    =generatenewmesh(Aloo,jacloo,slackloo,refl,tind,A,x(1:2*ncoils+NA(3)),ncoils,[stv(3) env(3)]/condi,[48 70]/condi);
    end
end

Jpp=sqrt((Apx(rangec,1:ncoils)*vmat*x(1:ncoils)).^2+...
         (Apy(rangec,1:ncoils)*vmat*x(1:ncoils)).^2+...
         (Apz(rangec,1:ncoils)*vmat*x(1:ncoils)).^2);
       max(Jpp)



%% cost function
jac=jac*10000;%volume at most 1
NA(3)=nnz(slack==1);
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
lb(stvar+1:stvar+ncoils)=0;
ub(stvar+1:stvar+ncoils)=Wmax;

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
Aeq = zeros(3,ncoils+NA(3)+ncoils);
Aeq(1,1:ncoils) = -Avector(1,1:ncoils,1)*condi;
Aeq=sparse(Aeq);
beq(1)=-50;
Aeq(2,stvar+1:end)=1;
beq(2)=Wmax;
%Aeq(3,1:ncoils) = -Avector(1,1:ncoils,2)*condi;
%beq(3)=-100;

Aeq(end+1,1:ncoils)=Apx(9848,1:ncoils)*vmat;
beq(end+1)=Jmax;
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
[x, fval, exitflag, output]  = intlinprog(f,intcon,sparse([Aine;Aine2]),sparse([bine;bine2]),Aeq,beq,lb,ub,xxx,options);
%%
if exitflag<=0
x=zeros([2*ncoils+NA(3) 1]);
break;
end
%%


end

initj=x(1:ncoils);
truequadratic=initj'*LL*initj;
zslack=sum(x(stvar+1:end-1))
initL1=max(abs(AL2*x(1:ncoils)))
pause(2)
end