ct=1; 
ss=[3.15 5.29 5.51 6.82 10.67 15.04 ];
ww=[183 114 107 103 55 40 86 102];
%for dept=[47 62 63 71]
    
for dept=[96 125 128 142 168 199]
load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si','peakJ','peakSR','Escalp','Escalp12','Ndep');
subplot(1,2,1),
plot(Wi,Escalp12)
hold on
xlabel('W (J)','Interpreter','latex');
ylabel('scalp $S_{1/2}$ $(cm^2)$','Interpreter','latex');

subplot(1,2,2),
plot(Wi,Escalp)
hold on
xlabel('W (J)','Interpreter','latex');
ylabel('peak E (V/m)','Interpreter','latex');

lentry{ct}=strcat('$ d_{1/2}$=',num2str(dept/10), '$(mm)$');
ct=ct+1;
end
legend(lentry{1},lentry{2},lentry{3},lentry{4},'Interpreter','latex');


%%

ct=1; 
ss=[3.15 5.29 5.51 6.82 10.67 15.04 ];
ww=[183 114 107 103 55 40 86 102];
for dept=[47 62 63 71 86 102]
load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si','peakJ','peakSR','Escalp','Escalp12','Ndep');
subplot(3,2,ct),
plot(Wi,Si)
hold on
xlabel('W (J)','Interpreter','latex');
ylabel('$ S_{1/\sqrt{2}}$ $(cm^2)$','Interpreter','latex');
title(strcat('$ d_{1/\sqrt{2}}$=',num2str(dept/10), '$(mm)$'),'Interpreter','latex');
scatter(ww(ct),ss(ct),'filled')
ct=ct+1;
end


ct=1;
for dept=[96 125 128 142  168 199]
load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si','peakJ','peakSR','Escalp','Escalp12','Ndep');
subplot(3,2,ct),
plot(Wi,Si)
hold on
xlabel('W (J)','Interpreter','latex');
ylabel('$ S_{1/\sqrt{2}}$ $(cm^2)$','Interpreter','latex');
ct=ct+1;
end