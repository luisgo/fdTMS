clear s
Ndep2=102;
load(strcat('diffdepthhat100results',num2str(Ndep2),'.mat'))
des=3;
for des=1:5
cur=100/sqrt(2)/Edecay{des}(1,1,Ndep2+1,1);
Edecay{des}=cur*Edecay{des};
Evol{des}=cur*Evol{des};

    for Ndep=1:400
S(Ndep)=10^-2*sum(Evol{des}(:)>=Edecay{des}(1,1,Ndep+1,1))/(Ndep/10);

end
subplot(2,1,1),
plot((0:399)/10,S)
xlabel('d (mm)')
xlim([0 10.2])
ylabel('S (cm^2)')
hold on
subplot(2,1,2),
plot((0:400)/10,squeeze(Edecay{des}(1,1,1:401,1)))
xlabel('d (mm)')
ylabel('V/m')
hold on

xlim([0 10.2])
end
legend('W=50','W=100','W=150','W=200','W=300')