dept=125;Wc=117;Sc=11.5;
subplot(3,2,1),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    plot(Wi,Si,'cyan','linewidth',2);
    hold on
    scatter(Wc,Sc,'red','filled');
    xlim([100 300]);
    title(strcat('d=',num2str(dept/10)));
    ylabel('S')
    subplot(3,2,2),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    plot(Wi,(1-Si/Sc)*100,'cyan','linewidth',2);
        ylabel('1-S/Sc (%)')
    xlim([100 300]);
    title(strcat('d=',num2str(dept/10)));
    
    dept=128;Wc=113;Sc=12.1;
subplot(3,2,3),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    plot(Wi,Si,'cyan','linewidth',2);
    hold on
    scatter(Wc,Sc,'red','filled');
    xlim([100 300]);   
    ylabel('S')
 
    title(strcat('d=',num2str(dept/10)));
    subplot(3,2,4),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    plot(Wi,(1-Si/Sc)*100,'cyan','linewidth',2);
    xlim([100 300]);
        ylabel('1-S/Sc (%)')
    title(strcat('d=',num2str(dept/10)));
    dept=142;Wc=107;Sc=15.2;
subplot(3,2,5),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    plot(Wi,Si,'cyan','linewidth',2);
    hold on
    scatter(Wc,Sc,'red','filled');
    xlim([100 300]);
       ylabel('S')
 
    title(strcat('d=',num2str(dept/10)));
    subplot(3,2,6),
    load(strcat('C:\Users\ljg24\Desktop\Coiloptcode\optimize\diffdepthhatresults',num2str(dept),'.mat'),'Wi','Si');
    plot(Wi,(1-Si/Sc)*100,'cyan','linewidth',2);
    xlim([100 300]);
    title(strcat('d=',num2str(dept/10)));
        ylabel('1-S/Sc (%)')
    