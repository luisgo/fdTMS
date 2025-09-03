clear all
syms y(t) x(t)
Dy=diff(y,t,1);
D2y=diff(y,t,2);
D3y=diff(y,t,3);
ode=D3y-Dy==exp(-3*t);
cond1 = y(0) == 23/24;
cond2 = Dy(0) == 3/24;
cond3 = D2y(0) == -9/24;
conds=[cond1 cond2 cond3];
dsolve(ode,conds)
ySol(t) = simplify(dsolve(ode,conds))
