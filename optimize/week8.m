h=1/2;
Nsteps=round(40/h);
R=5;C=1;v0=4;
vin=@(t) t.*(t>=0);
Dx= @(t,vt) (vin(t)-vt)/(R*C);
t=zeros([Nsteps 1]);
v=zeros([Nsteps 1]);
dv=zeros([Nsteps 1]);
v(1)=v0;
for i=1:Nsteps-1
    t(i)=(i-1)*h;
    dv(i)=Dx(t(i),v(i));%diffeq
    v(i+1)=v(i)+h*dv(i);%approximation
end
 t(Nsteps)=(Nsteps-1)*h;
 
plot(t,vin(t),'linewidth',2);
hold on
plot(t,v,'linewidth',2);
plot(t,dv,'linewidth',2);
legend('v_{in}','v_c(t)','dV_c(t)/dt')