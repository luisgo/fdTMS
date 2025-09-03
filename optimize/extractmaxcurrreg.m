%%
clear xv yv
harmfile='../newwireharmonics.mat';
dept=101;Wc=107;Sc=15.2;
Ncc=24;
coilid=1;
Ncc=24;
Nco=5;
dept=[101]
    dept2=[49 49 49 49 49 49 49 49];
filename='optimT3no100';
filenamesav=strcat(filename,'d',num2str(dept),'ncoils',num2str(Ncc));
designfile=strcat(filenamesav,'lmax30.mat');
for coilid=1:Nco;

load(harmfile);
load(designfile);
load  indevals.mat inde Jparty ASR;
xv(coilid)=max(Jparty(:,1:Ncc)*Jcalc2{coilid})
yv(coilid)=max(ASR(:,1:Ncc)*Jcalc2{coilid})
end
%%
mag=(yv./xv/5);
%%
subplot(3,1,1),plot(1:Nco,(yv./xv/5))
xlabel('maximum derivative')
ylabel('gap')

subplot(3,1,2),plot(1:Nco,Si)
xlabel('maximum derivative')
ylabel('Si')

subplot(3,1,3),plot(1:Nco,Wi)
xlabel('maximum derivative')
ylabel('Wi')
%%
load(harmfile);
load(designfile);
coilid=3;
Att=Apx(:,1:Ncc)*Jcalc2{coilid};
Att=Att(:,:);
Att=Asym*Atrans(:,1:Ncc)*Jcalc2{coilid};
trisurf(t2psym,psym(:,1),psym(:,2),psym(:,3),Att(:),'edgealpha',0);

inde=1:numel(Att);
lo=(loc(:,1)<0).*(loc(:,2)<0).*(loc(:,2)>=-.018).*(loc(:,1)>=-.02);
inde=inde(lo>0);
hold on
trisurf(t2psym(inde,:),psym(:,1),psym(:,2),psym(:,3),'facealpha',.3,'facecolor','red','edgealpha',0);

col=zeros([3*numel(inde) 1]);
row=zeros([3*numel(inde) 1]);
val=zeros([3*numel(inde) 1]);
for i=1:numel(inde)
col(3*i-2:3*i)=i;
ps=t2psym(inde(i),:);
row(3*i-2:3*i)=ps(:);
%shift origin to center of triangle
% pv=psym(ps,:);
% pcen=(pv(1,:)+pv(2,:)+pv(3,:))/3;
% pv(:,1)=pv(:,1)-pcen(1);
% pv(:,2)=pv(:,2)-pcen(2);
% pv(:,3)=pv(:,3)-pcen(3);
% mat=pinv(pv,10^-8*(norm(pv)));
% val(3*i-2:3*i)=(mat(2,:)*(eye(3)-ones([3 3])/3))';
mat=inv(psym(ps,:));
val(3*i-2:3*i)=(mat(2,:))';
end
Jparty=sparse(col,row,val,numel(inde),numel(psym)/3)*Asym*Atrans(:,:);
Jtt=Jparty(:,1:Ncc)*Jcalc2{coilid};
trisurf(t2psym(inde,:),psym(:,1),psym(:,2),psym(:,3),Jtt,'edgealpha',0);
ASR=Asym*Atrans(:,1:Ncc);




save indevals.mat inde Jparty ASR;