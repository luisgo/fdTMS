ct=1;
for dept=[96 128 142 200]
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhat100results',num2str(dept),'.mat'));
subplot(2,2,ct),
plot(Wi,Escalp);
xlabel(strcat('depth=',num2str(dept)))
ct=ct+1
end