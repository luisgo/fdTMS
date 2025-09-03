clear all;
harmfile='../newwireharmonics.mat';
load 'C:\Users\ljg24\Desktop\Coiloptcode\generatecoilmesh\firstordercorr_middle.mat'
concentric=10;
[designfile resultsfile Ncc]=getfilename('thick100',4,3,coilid);
load(designfile)
if numel(Jcalc2)==1
coilid=1;
end
load(resultsfile,'Si','Wi');

[xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,coilid,concentric(ct2),Ncc);
fb80rs=[];fb80js=[];
clen=0;
resolu=.0001;
ct=1;
wid=4.067098615315901e-03/2;
R=0;
for i=1:numel(xcoil)
    clear nhat;
    nhat(:,1)=nx(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
    nhat(:,2)=ny(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
    nhat(:,3)=nz(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
 

ra=[xcoil{i}([1:end ])'+wid*nhat(:,1),...
ycoil{i}([1:end ])'+wid*nhat(:,2),...
zcoil{i}([1:end ])'+wid*nhat(:,3)];
rb=[xcoil{i}([1:end ])'-wid*nhat(:,1),...
ycoil{i}([1:end ])'-wid*nhat(:,2),...
zcoil{i}([1:end ])'-wid*nhat(:,3)];      
len=cumsum(sqrt([0;sum((ra(2:end,:)-ra(1:end-1,:)).^2,2)]));
R=R+len(end);
xa = interp1(len,ra(:,1),0:.0002:len(end));
ya = interp1(len,ra(:,2),0:.0002:len(end));
za = interp1(len,ra(:,3),0:.0002:len(end));
ra=[xa(:),ya(:),za(:)];

len=cumsum(sqrt([0;sum((rb(2:end,:)-rb(1:end-1,:)).^2,2)]));
R=R+len(end);
xb = interp1(len,rb(:,1),0:.0002:len(end));
yb = interp1(len,rb(:,2),0:.0002:len(end));
zb = interp1(len,rb(:,3),0:.0002:len(end));
rb=[xb(:),yb(:),zb(:)];


fb80rs=cat(1,fb80rs,...
    (ra(2:end,:)+ra(1:end-1,:))/2,...
    (rb(2:end,:)+rb(1:end-1,:))/2);       
fb80js=cat(1,fb80js,...
    (ra(2:end,:)-ra(1:end-1,:)),...
    (rb(2:end,:)-rb(1:end-1,:)));
     ct=ct+1;
end
R80=R;

rv=fb80rs';Jdc=fb80js';

%%
concentric=8;

    close all
dept=131;Wc=107;Sc=15.2;
Ncc=24;
coilid=4;
designfile=strcat('optimT3no100d',num2str(dept),'ncoils24lmax30.mat');
    load(strcat('.\diffdepthhatresults',num2str(dept),'.mat'),'Si','Wi');
    ct=1;
i=4;

[xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,coilid,concentric,Ncc);
fb65rs=[];fb65js=[];
clen=0;
resolu=.0001;
ct=1;
R=0;
for i=1:numel(xcoil)
    
    clear nhat;
    nhat(:,1)=nx(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
    nhat(:,2)=ny(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
    nhat(:,3)=nz(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
 

ra=[xcoil{i}([1:end ])'+wid*nhat(:,1),...
ycoil{i}([1:end ])'+wid*nhat(:,2),...
zcoil{i}([1:end ])'+wid*nhat(:,3)];
rb=[xcoil{i}([1:end ])'-wid*nhat(:,1),...
ycoil{i}([1:end ])'-wid*nhat(:,2),...
zcoil{i}([1:end ])'-wid*nhat(:,3)];      
len=cumsum(sqrt([0;sum((ra(2:end,:)-ra(1:end-1,:)).^2,2)]));
xa = interp1(len,ra(:,1),0:.0002:len(end));
ya = interp1(len,ra(:,2),0:.0002:len(end));
za = interp1(len,ra(:,3),0:.0002:len(end));
ra=[xa(:),ya(:),za(:)];
R=R+len(end);
len=cumsum(sqrt([0;sum((rb(2:end,:)-rb(1:end-1,:)).^2,2)]));
xb = interp1(len,rb(:,1),0:.0002:len(end));
yb = interp1(len,rb(:,2),0:.0002:len(end));
zb = interp1(len,rb(:,3),0:.0002:len(end));
rb=[xb(:),yb(:),zb(:)];

R=R+len(end);
fb65rs=cat(1,fb65rs,...
    (ra(2:end,:)+ra(1:end-1,:))/2,...
    (rb(2:end,:)+rb(1:end-1,:))/2);       
fb65js=cat(1,fb65js,...
    (ra(2:end,:)-ra(1:end-1,:)),...
    (rb(2:end,:)-rb(1:end-1,:)));
     ct=ct+1;
end

R65=R;
%%
R=0;

concentric=8;

    close all
dept=101;Wc=107;Sc=15.2;
Ncc=24;
coilid=4;
designfile=strcat('optimT3no100d',num2str(dept),'ncoils24lmax30.mat');
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Si','Wi');
ct=1;
i=5;


[xcoil,ycoil,zcoil,zzcoil,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,coilid,concentric,Ncc);
fb35rs=[];fb35js=[];
clen=0;
resolu=.0001;
ct=1;
for i=1:numel(xcoil)
    
    clear nhat;
    nhat(:,1)=nx(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
    nhat(:,2)=ny(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
    nhat(:,3)=nz(xcoil{i}([1:end ]),ycoil{i}([1:end ]));
 

ra=[xcoil{i}([1:end ])'+wid*nhat(:,1),...
ycoil{i}([1:end ])'+wid*nhat(:,2),...
zcoil{i}([1:end ])'+wid*nhat(:,3)];
rb=[xcoil{i}([1:end ])'-wid*nhat(:,1),...
ycoil{i}([1:end ])'-wid*nhat(:,2),...
zcoil{i}([1:end ])'-wid*nhat(:,3)];      
len=cumsum(sqrt([0;sum((ra(2:end,:)-ra(1:end-1,:)).^2,2)]));
xa = interp1(len,ra(:,1),0:.0002:len(end));
ya = interp1(len,ra(:,2),0:.0002:len(end));
za = interp1(len,ra(:,3),0:.0002:len(end));
ra=[xa(:),ya(:),za(:)];
R=R+len(end);
len=cumsum(sqrt([0;sum((rb(2:end,:)-rb(1:end-1,:)).^2,2)]));
xb = interp1(len,rb(:,1),0:.0002:len(end));
yb = interp1(len,rb(:,2),0:.0002:len(end));
zb = interp1(len,rb(:,3),0:.0002:len(end));
rb=[xb(:),yb(:),zb(:)];
R=R+len(end);

fb35rs=cat(1,fb35rs,...
    (ra(2:end,:)+ra(1:end-1,:))/2,...
    (rb(2:end,:)+rb(1:end-1,:))/2);       
fb35js=cat(1,fb35js,...
    (ra(2:end,:)-ra(1:end-1,:)),...
    (rb(2:end,:)-rb(1:end-1,:)));
     ct=ct+1;
end
R35=R;
addpath('D:\Coiloptcode\generatecoilmesh\figure8coilmodelgenerator\wirecode');


save twolayerdes.mat fb35rs fb35js fb65rs fb65js fb80rs fb80js R35 R65 R80;