function Ahat=generateL1mat_adapt(that1,that2,Ax,Ay,Az,rangec,Nmodes,Nang,pts)
Ahat=zeros([Nang*numel(pts),Nmodes]);
ct=1;
for ipp=1:numel(pts)
    ipt=pts(ipp);
for j=1:Nang
    th=2*pi/Nang*(j-1/2);
for i=1:Nmodes
Ahat(ct,i)=...
(that1(ipt,1).*Ax(ipt,i)+that1(ipt,2).*Ay(ipt,i)+that1(ipt,3).*Az(ipt,i))*sin(th)+...;
(that2(ipt,1).*Ax(ipt,i)+that2(ipt,2).*Ay(ipt,i)+that2(ipt,3).*Az(ipt,i))*cos(th);
end
ct=ct+1;
end
end