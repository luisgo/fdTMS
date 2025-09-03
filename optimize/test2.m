clear all
p1=@(t) (t>=-0.5).*(t<=0.5);
p2=@(t) (t>=-0.3).*(t<=0.3);
ct=1;tau=-1:.01:1;tau2=-1:.001:1;
for t=-1:.2:1
    int(ct)=(sum(p1(tau2).*p2(t-tau2)))*.001;
    subplot(4,3,ct)
    plot(tau,p1(tau),'red');
    hold on
    plot(tau,p2(t-tau),'yellow');
    xlabel('\tau')
    title(strcat('t=',num2str(t)))
    scatter(t,int(ct),'filled','blue')
        ct=ct+1;

end
    subplot(4,3,ct),
    plot(-1:.2:1,int,'o-');ylim([0 1]);
    hold on
   scatter(-1:.2:1,int,'filled','blue')
xlabel('t')