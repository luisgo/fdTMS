%% problem 2a
clear all
kn=4*10^-3; VTh=2;
syms ID;
VDS=-375*ID+10;
VGS=-125*ID+8;
IDsaturation=kn/2*(VGS-VTh)^2;
c=sym2poly(IDsaturation-ID);
%IDtriode=kn*(VGS-VTh)*VDS-kn/2*VDS^2;
%c=sym2poly(IDtriode-ID)
ID=roots(c);
VDS=-375*ID+10;
VGS=-125*ID+8;
ID(VGS>VTh)
VDS(VGS>VTh)
VGS(VGS>VTh)
%% problem 2b
clear all
kn=4*10^-3; VTh=-2;
syms ID;
VDS=-3000*ID-10;
VGS=-1000*ID-6.5;
%IDsaturation=-kn/2*(VGS-VTh)^2;
%c=sym2poly(IDsaturation-ID)
IDtriode=-kn*(VGS-VTh)*VDS+kn/2*VDS^2;
c=sym2poly(IDtriode-ID);
ID=roots(c);
VDS=-3000*ID-10;
VGS=-1000*ID-6.5;
ID(VDS>VGS-VTh)
VDS(VDS>VGS-VTh)
VGS(VDS>VGS-VTh)
%% problem 3a
Rss=500;
R2=2000;
ID=.003;
VGS=-2*sqrt(3)-1
VDS=VGS+1
VS=3-ID*Rss
VG=VS+VGS
R1=(6*R2)/(VG+3)-R2
RD=-(-VDS-6+.003*Rss)/.003

%% problem 3b
Rss=500;
R2=2000;
ID=.003;
VGS=2*sqrt(3)-1
VS=ID*Rss
VG=VS+VGS
R1=(10*R2)/VG-R2
VDS=VGS+1;
RD=(10-VDS-.003*Rss)/.003
%HW 1 soln
%RD=1500;
%VDS=10-ID*(Rss+RD);
%% problem 4
VDS=1;VGS=5;R=1000;VTh=1;
%triode region
ID=(5-VDS)/R;
kn=ID/((VGS-VTh)*VDS-1/2*VDS^2)

%% problem 4b
Vp=0:.001:10;
VGS=0; VTh=-2;kn=1;
I=kn*(VGS-VTh)*Vp-1/2*Vp.^2
I(Vp>=VGS-VTh)=kn*(VGS-VTh)^2/2;
plot(Vp,I)
xlabel('V_p (V)');
ylabel('I (mA)');

%% problem 4c
Vp=0:.001:10;
VGS=0; VTh=-2;kn=1;
I=kn*(VGS-VTh)*Vp-1/2*Vp.^2
I(Vp>=VGS-VTh)=kn*(VGS-VTh)^2/2;
plot(Vp,I)
xlabel('V_p (V)');
ylabel('I (mA)');
%% problem 5b


%% problem 2b

clear all
kn=4*10^-3; VTh=-2;
syms ID;
VDS=-3000*ID-10;
VGS=-1000*ID-6.5;
%IDsaturation=-kn/2*(VGS-VTh)^2;
%c=sym2poly(IDsaturation-ID)
IDtriode=-kn*(VGS-VTh)*VDS+kn/2*VDS^2;
c=sym2poly(IDtriode-ID);
ID=roots(c);
VDS=-3000*ID-10;
VGS=-1000*ID-6.5;
ID(VDS>VGS-VTh)
VDS(VDS>VGS-VTh)
VGS(VDS>VGS-VTh)
