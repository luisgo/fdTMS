function [no2,ou]=pt2line(a,c1,c2);
nhat=[-(c2(2)-c1(2)) (c2(1)-c1(1))];
no=norm(nhat);
nhat=nhat/no;
no2=sum((c1-a).*nhat);
ou=a+no2*nhat;