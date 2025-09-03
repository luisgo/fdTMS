clear all
figure
 yv=[8.83;14.8;22.85];
xv=1:10:500;
xv3=300:10:800;
xv2=[206;105;57]
dv=[10;14;17;20];
1-[9.044 8 7.668]/14.1
enval=3;
SP=.1
SV=.1
ppp=0;
pad=.1
mar=.1
for i=1:enval
if ppp==1;
    subaxis(1,3,i,'SpacingHoriz',  SP, 'SpacingVert', SV, 'padding', pad, 'margin', mar);
    else
        subplot(1,3,i),
end
 hold on
% title(strcat('d=',num2str(dv(i)),'mm'))
     scatter(xv2(i),yv(i),'r','filled');
     if i==1
     text(xv2(i),yv(i),' \leftarrow 25 mm Figure-8','FontName','Times','FontSize',20);
     elseif i==2
     text(xv2(i),yv(i),' \leftarrow 70 mm Figure-8','FontName','Times','FontSize',20);
     elseif i==3
     text(xv2(i),yv(i),' \leftarrow Double cone','FontName','Times','FontSize',20);
     end
     end
% for i=2:enval
% subplot(2,2,i),
% title(strcat('d=',num2str(dv(i)),'mm'))
%     scatter(xv(i),yv(i),'r','filled');
% end
    hold on
    col{1}='X-b';
    col{2}='X-g';
    col{3}='X-black';
 for kkk=1:3
 for i=1:enval
     if kkk==1
     load('C:\Users\ljg24\Desktop\optimisationcodes\spherecoildesfiles\energ.mat','Wi','Si')
     elseif kkk==2
     load('C:\Users\ljg24\Desktop\optimisationcodes\hemispherecoildesfiles\energ.mat','Wi','Si')
     elseif kkk==3
     load('C:\Users\ljg24\Desktop\optimisationcodes\planarcoildesfiles\energ.mat','Wi','Si')   
     end
 Si{i}=Si{i}(isfinite(Wi{i}));
 Wi{i}=Wi{i}(isfinite(Wi{i}));
 Wi{i}=Wi{i}(Si{i}~=0);
 Si{i}=Si{i}(Si{i}~=0);
 [~,x]=sort(Wi{i});
 Si{i}=Si{i}(x);
 Wi{i}=Wi{i}(x);    

if ppp==1;
    subaxis(2,2,i,'SpacingHoriz',  SP, 'SpacingVert', SV, 'padding', pad, 'margin', mar);
    else
        subplot(1,3,i),
end
    k = convhull(Wi{i}(1:end),Si{i}(1:end));
    %k=1:numel(Wi{i});
  %  if i==4 && kkk==2
  % k=k(1:end-3);
  k=unique(k,'stable');
  %  else
   % k=k(1:end);
   % end
    hold on
    %trim bow
    while Wi{i}(k(end))<Wi{i}(k(end-1))
        k=k(1:end-1);
    end
    plot(Wi{i}(k),Si{i}(k),strcat(col{kkk}),'linewidth',2)
    k=1:numel(Wi{i});

 %  scatter(Wi{i}(k),Si{i}(k),col{kkk}(3),'filled')
%plot(Wi{i}(k),Si{i}(k),'X-r','linewidth',2)
xlabel('Energy (J)')
ylabel('Spread (cm^2)')
set(gca,'FontName','Times','FontSize',24)
xlim([10 500])

%ylim([floor(1.21*min(Si{i}(k))) ceil(max(Si{i}(k)))])
 end
 end

if ppp==1;
    subaxis(2,2,1,'SpacingHoriz',  SP, 'SpacingVert', SV, 'padding', pad, 'margin', mar);
    else
        subplot(1,3,1),
end
    ylim([3 10])
    set(gca,'ytick',[3 6 8 10])
    grid on

if ppp==1;
    subaxis(2,2,2,'SpacingHoriz',  SP, 'SpacingVert', SV, 'padding', pad, 'margin', mar);
    else
        subplot(1,3,2),
end
ylim([7 16])
    grid on

if ppp==1;
    subaxis(2,2,3,'SpacingHoriz',  SP, 'SpacingVert', SV, 'padding', pad, 'margin', mar);
    else
        subplot(1,3,3),
end 
    ylim([12 24])
    set(gca,'ytick',[12 16 20 24])
    grid on
dept=100;
subplot(1,3,1),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
plot(Wi,Si,'cyan','linewidth',2);
dept=140;
subplot(1,3,2),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
plot(Wi,Si,'cyan','linewidth',2);
dept=170;

    subplot(1,3,3),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    plot(Wi,Si,'cyan','linewidth',2);
     legend('50-Coils','Sphere','Half-sphere','Square','Hat')

 %%
 figure
 j=enval;
 if j==1
     
kv{1}=[1 9 14 16];
load 'C:\qTMS Coil - PCB\planenogaparebranchandbound11front20ncoils100.mat'
Jcalc3=Jcalc2;
load 'C:\qTMS Coil - PCB\planenogaparelowpower11front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
load 'C:\qTMS Coil - PCB\planenogaparebranchendbound11front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
load 'C:\qTMS Coil - PCB\planenogaparerunfast11front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
Jcalc2=Jcalc3;
 elseif j==2
kv{2}=[1  7 10 12];
load 'C:\qTMS Coil - PCB\planenogaparebranchandbound14front20ncoils100.mat'
Jcalc3=Jcalc2;
load 'C:\qTMS Coil - PCB\planenogaparelowpower14front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
load 'C:\qTMS Coil - PCB\planenogaparebranchendbound14front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
load 'C:\qTMS Coil - PCB\planenogaparerunfast14front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
Jcalc2=Jcalc3;
 elseif j==3
kv{3}=[ 1  5 6 7];
load 'C:\qTMS Coil - PCB\planenogaparebranchandbound17front20ncoils100.mat'
Jcalc3=Jcalc2;
load 'C:\qTMS Coil - PCB\planenogaparelowpower17front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
load 'C:\qTMS Coil - PCB\planenogaparebranchendbound17front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
load 'C:\qTMS Coil - PCB\planenogaparerunfast17front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
Jcalc2=Jcalc3;
 elseif j==4
kv{4}=[1  4 7 9];
load 'C:\qTMS Coil - PCB\planenogaparebranchandbound20front20ncoils100.mat'
Jcalc3=Jcalc2;
load 'C:\qTMS Coil - PCB\planenogaparelowpower20front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
load 'C:\qTMS Coil - PCB\planenogaparebranchendbound20front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
load 'C:\qTMS Coil - PCB\planenogaparerunfast20front20ncoils100.mat'
Jcalc3=cat(2,Jcalc2,Jcalc3);
Jcalc2=Jcalc3;
 end
 clear Jcalc3;
 ct=1;
 Wi2=Wi;Si2=Si;
 load results1.mat Wi Si  
 for i=1:numel(Jcalc2)
     if isfinite(Wi{enval}(i))==1;
 Jcalc3{ct}=Jcalc2{i};
 ct=ct+1;
     end
 end
 Wi=Wi2;Si=Si2;
 Jcalc2=Jcalc3;
   
 for i=1:4
 kval(i)=k(kv{enval}(i));
  ncoils=100
  Jcalc=Jcalc2{x(k(kv{j}(i)))}(1:ncoils);
  Jv=Jcalc(:)/(2*pi*3000);
  [winding2{i},Idrive{i},col{i}]=drawcoilplane(eye(ncoils),Jcalc,sqrt(ncoils));
   end
  
   close all
for i=1:4
  winding=winding2{i};
  Nc=numel(winding(1,1,:));
  js=zeros([3 300*Nc]);
  rs=zeros([3 300*Nc]);
for j=1:Nc
js(:,300*(j-1)+1:300*j)=winding(:,[2:end,1],j)-winding(:,1:end,j);
rs(:,300*(j-1)+1:300*j)=(winding(:,[2:end,1],j)+winding(:,1:end,j))/2;
end
load 'C:\qTMS Coil - PCB\interppoints.mat' tri p
p=.069*p;
rhos=zeros(numel(rs(1,:)),1);
ks=zeros(size(rs));
Efield=computefields(rs',js',ks',rhos,2*pi*3000,15,.086,p);
Efield=sqrt(sum(abs(Efield).^2,2));
subplot(2,2,i),
trisurf(tri,p(:,1),p(:,2),p(:,3),Efield,'Facecolor','Interp','edgealpha',0)
hold on
for j=1:Nc
plot3(rs(1,300*(j-1)+1:300*j),rs(2,300*(j-1)+1:300*j),rs(3,300*(j-1)+1:300*j),col{i}(j),'linewidth',1)
end
view([0 0 1])
  title(strcat('(',num2str(round(Wi{enval}(kval(i)))),',',num2str(round(10*Si{enval}(kval(i)))/10),')'))
axis equal
axis off
  set(gca,'FontName','Times','FontSize',24)

end