%%double cone
subplot(2,3,1),
plot(Wi,Ssq,'o-','linewidth',3)
hold on; scatter(33,8.8,'filled')
xlabel('energy (J)');
ylabel('S_{1/sqrt(2)}')
title('depth=7.9 mm')
ylim([4 9])
subplot(2,3,2),
plot(Wi,Ssq/8.8,'o-','linewidth',3)
hold on; scatter(33,1,'filled')
xlabel('energy (J)');
ylabel('S_{1/sqrt(2)}')
title('depth=7.9 mm')
ylim([0.5 1]);
subplot(2,3,3),
plot(Wi,Ndepsq/10,'o-','linewidth',3)
hold on;
xlabel('energy (J)');
ylabel('D_{1/sqrt(2)}')
title('depth=7.9 mm')
hold on; scatter(33,7.9,'filled')


subplot(2,3,4),
plot(Wi,Si,'o-','linewidth',3)
hold on; scatter(33,19.8,'filled')
xlabel('energy (J)');
ylabel('S_{1/2}')
title('depth=15.7 mm')
ylim([12 20]);
subplot(2,3,5),
plot(Wi,Si/19.8,'o-','linewidth',3)
hold on; scatter(33,1,'filled')
xlabel('energy (J)');
ylabel('S_{1/2}')
title('depth=15.7 mm')
ylim([0.6 1]);
%% fig8

subplot(2,3,1),
plot(Wi,Ssq,'o-','linewidth',3)
hold on; scatter(125,5.8,'filled')
xlabel('energy (J)');
ylabel('S_{1/sqrt(2)}')
title('depth=6.5 mm')
ylim([2.8 6])
subplot(2,3,2),
plot(Wi,Ssq/5.8,'o-','linewidth',3)
hold on; scatter(125,1,'filled')
xlabel('energy (J)');
ylabel('S_{1/sqrt(2)}')
title('depth=6.5 mm')
ylabel([0.5 1]);
subplot(2,3,3),
plot(Wi,Ndepsq/10,'o-','linewidth',3)
xlabel('energy (J)');
ylabel('D_{1/sqrt(2)}')
title('depth=6.5 mm')
hold on
hold on; scatter(125,6.5,'filled')
subplot(2,3,4),
plot(Wi,Si,'o-','linewidth',3)
hold on; scatter(125,12.8,'filled')
xlabel('energy (J)');
ylabel('S_{1/2}')
title('depth=13.1 mm')
ylim([8 13]);
subplot(2,3,5),
plot(Wi,Si/12.8,'o-','linewidth',3)
hold on; scatter(125,1,'filled')
xlabel('energy (J)');
ylabel('S_{1/2}')
title('depth=13.1 mm')
ylim([.6 1]);

%% small

subplot(2,3,1),
plot(Wi,Ssq,'o-','linewidth',3)
hold on; scatter(215,3.37,'filled')
xlabel('energy (J)');
ylabel('S_{1/sqrt(2)}')
title('depth=4.9 mm')
subplot(2,3,2),
plot(Wi,Ssq/3.37,'o-','linewidth',3)
hold on; scatter(215,1,'filled')
xlabel('energy (J)');
ylabel('S_{1/sqrt(2)}')
title('depth==4.9 mm')
subplot(2,3,3),
plot(Wi,Ndepsq,'o-','linewidth',3)
xlabel('energy (J)');
ylabel('D_{1/sqrt(2)}')
title('depth=4.9 mm')

hold on; scatter(215,4.9,'filled')
subplot(2,3,4),
plot(Wi,Si,'o-','linewidth',3)
hold on; scatter(215,7.5,'filled')
xlabel('energy (J)');
ylabel('S_{1/2}')
title('depth=10.1 mm')
subplot(2,3,5),
plot(Wi,Si/7.5,'o-','linewidth',3)
hold on; scatter(230,1,'filled')
xlabel('energy (J)');
ylabel('S_{1/2}')
title('depth=10.1 mm')