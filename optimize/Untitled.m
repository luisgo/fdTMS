figure
xlabel('Re {S}','Interpreter','latex','fontsize',20)
ylabel('Im {S}','Interpreter','latex','fontsize',20)
xlim([-1 1]);
ylim([-1 1]);

xL = xlim;
yL = ylim;
line([0 0], yL,'linewidth',2);  %x-axis
line(xL, [0 0],'linewidth',2);  %y-axis
set(gca,'fontsize',20,'fontname','tex')