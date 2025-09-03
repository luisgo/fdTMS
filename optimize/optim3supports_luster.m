function optim3supports_cluster(iicluster)
if iicluster<=3
harmfile='../halfsphereharm2.mat';
Ncc=22;
ct=1;
 dept=[101 131 157]
dept=dept(iicluster);
filename='paperover2optimshalfsphere2';
filenamesav=strcat(filename,'d',num2str(dept),'ncoils',num2str(Ncc));
designfile=strcat(filenamesav,'lmax30.mat');
ener=[50 100 150 200 250 300];
peakSR2=10^9*ones(size(ener));
%Jmax=10^6*[1.5832    1.9988    3.2843    3.2195    4.0595]*.9;%125
%Jmax=10^6*[3.2052    5.1185    5.6964    7.1927    9.3359]*.9;%96
Jmax=10^12*ones(size(ener));%142
gap=10^30*ones(size(ener));
 %ncoils = runoptiminogaL12layer2(dept,ener(:),Jmax(:),gap(:),Ncc,harmfile,filename);
%ncoils = runoptiminogatwo([dept,dept2(ct)],ener,Jmax,peakSR,Ncc,harmfile,filename);
 ncoils = runoptiminoga2(dept,ener,Jmax,peakSR2,Ncc,harmfile,filename);
 %ncoils = runoptiminoga(dept,ener,Ncc,harmfile,filename);
close all
addpath('C:\Users\ljg24\Desktop\postproc\determinetradeoff\')

addpath('..\determinetradeoff');

[Wi,Si,peakJ,peakSR,Escalp,Escalp12,Ndep,Evol,Edecay,Escalp,Ssq,Ndepsq,Sd,Wd]=anycoildeterminetradeoff_f(harmfile,designfile,dept,Ncc);
save(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\PAPERhalfsphere2WOUT1002',num2str(dept),'.mat'),'Wi','Si','Ssq','peakJ','peakSR','Escalp','Escalp12','Ndep','Evol','Edecay','Escalp','Ndepsq','Sd','Wd');
peakSR
ct=ct+1;
elseif iicluster<=6
harmfile='../spheremharm.mat';
Ncc=64;
%for dept=[49 65 79]
ct=1;
dept=[101 131 157]
dept=dept(iicluster-3);
filename='paperover2optimssphere';
filenamesav=strcat(filename,'d',num2str(dept),'ncoils',num2str(Ncc));
designfile=strcat(filenamesav,'lmax30.mat');
ener=[50 100 150 200 250 300];
peakSR2=10^9*ones(size(ener));
%Jmax=10^6*[1.5832    1.9988    3.2843    3.2195    4.0595]*.9;%125
%Jmax=10^6*[3.2052    5.1185    5.6964    7.1927    9.3359]*.9;%96
Jmax=10^12*ones(size(ener));%142
gap=10^30*ones(size(ener));
 %ncoils = runoptiminogaL12layer2(dept,ener(:),Jmax(:),gap(:),Ncc,harmfile,filename);
%ncoils = runoptiminogatwo([dept,dept2(ct)],ener,Jmax,peakSR,Ncc,harmfile,filename);
 ncoils = runoptiminoga2(dept,ener,Jmax,peakSR2,Ncc,harmfile,filename);
 %ncoils = runoptiminoga(dept,ener,Ncc,harmfile,filename);
close all
addpath('C:\Users\ljg24\Desktop\postproc\determinetradeoff\')

addpath('..\determinetradeoff');

[Wi,Si,peakJ,peakSR,Escalp,Escalp12,Ndep,Evol,Edecay,Escalp,Ssq,Ndepsq,Sd,Wd]=anycoildeterminetradeoff_f(harmfile,designfile,dept,Ncc);
save(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\PAPERsphereWOUT1002',num2str(dept),'.mat'),'Wi','Si','Ssq','peakJ','peakSR','Escalp','Escalp12','Ndep','Evol','Edecay','Escalp','Ndepsq','Sd','Wd');
peakSR
ct=ct+1;
elseif iicluster<=9
%%
harmfile='../squareharm2.mat';
Ncc=40;
%for dept=[49 65 79]
ct=1;
dept=[101 131 157]
dept=dept(iicluster-6);
filename='paperover2optimshalfsquare';
filenamesav=strcat(filename,'d',num2str(dept),'ncoils',num2str(Ncc));
designfile=strcat(filenamesav,'lmax30.mat');
ener=[50 100 150 200 250 300];
peakSR2=10^9*ones(size(ener));
%Jmax=10^6*[1.5832    1.9988    3.2843    3.2195    4.0595]*.9;%125
%Jmax=10^6*[3.2052    5.1185    5.6964    7.1927    9.3359]*.9;%96
Jmax=10^12*ones(size(ener));%142
gap=10^30*ones(size(ener));
 %ncoils = runoptiminogaL12layer2(dept,ener(:),Jmax(:),gap(:),Ncc,harmfile,filename);
%ncoils = runoptiminogatwo([dept,dept2(ct)],ener,Jmax,peakSR,Ncc,harmfile,filename);
 ncoils = runoptiminoga2(dept,ener,Jmax,peakSR2,Ncc,harmfile,filename);
 %ncoils = runoptiminoga(dept,ener,Ncc,harmfile,filename);
close all
addpath('C:\Users\ljg24\Desktop\postproc\determinetradeoff\')

addpath('..\determinetradeoff');

[Wi,Si,peakJ,peakSR,Escalp,Escalp12,Ndep,Evol,Edecay,Escalp,Ssq,Ndepsq,'Sd','Wd']=anycoildeterminetradeoff_f(harmfile,designfile,dept,Ncc);
save(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\PAPERhalfsquareWOUT1002',num2str(dept),'.mat'),'Wi','Si','Ssq','peakJ','peakSR','Escalp','Escalp12','Ndep','Evol','Edecay','Escalp','Ndepsq','Sd','Wd');
peakSR
ct=ct+1;
end

end