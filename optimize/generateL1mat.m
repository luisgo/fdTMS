function Ahat=generateL1mat(that1,that2,Ax,Ay,Az,rangec,Nmodes,Nang)
Ahat=zeros([Nang*numel(that1(:,1)),Nmodes]);
for j=1:Nang
    ran=(j-1)*numel(that1(:,1))+1:j*numel(that1(:,1));
    th=2*pi/Nang*(j-1/2);
for i=1:Nmodes
Ahat(ran,i)=...
(that1(:,1).*Ax(rangec,i)+that1(:,2).*Ay(rangec,i)+that1(:,3).*Az(rangec,i))*sin(th)+...;
(that2(:,1).*Ax(rangec,i)+that2(:,2).*Ay(rangec,i)+that2(:,3).*Az(rangec,i))*cos(th);
end
end