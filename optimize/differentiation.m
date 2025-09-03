syms x xp y yp

r(x,y,xp,yp)=((x-xp)^2+(y-yp)^2)^(1/2);
f(x,y,xp,yp)=(x-xp)*log(r(x,y,xp,yp))+...
    (y-yp)*atan((x-xp)/(y-yp))...
    -(x-xp);
diff(f,xp)
simplify(diff(f,xp))

simplify(diff(f,yp))
%%
p=rand([10 2]);
k=convhull(p);
scatter(p(k,1),p(k,2),'filled')
hold on
plot(p(k,1),p(k,2))
axis off
for i=1:numel(k)-1
    text(p(k(i),1),p(k(i),2)+.05,num2str(i))
end

