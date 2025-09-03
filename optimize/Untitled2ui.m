subplot(2,2,1),
scatter([-2 -2 -1],[1 -1 0],'x')
hold on
scatter([0],[0],'o')
subplot(2,2,2),
scatter([2 2 1],[1, -1 0],'x')
hold on
scatter([0],[0],'o')

subplot(2,2,3),
scatter([-1 -2],[0 0],'x')
hold on
scatter([2],[2],'o')

subplot(2,2,4),
scatter([-1 2],[0 0],'x')
hold on
scatter([0 2],[0 2],'o')
for i=1:4
subplot(2,2,i),
xlim([-5 5]);
ylim([-5 5]);
grid on
end