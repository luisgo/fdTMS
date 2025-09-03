function c=costfnall(clevs,Eideal,pmod2,t2pmod2,Attmod2,Jval,Nv,nhh)

levs=cat(1,clevs(1:Nv(1)),-clevs(1:Nv(1)));
[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{1}-.002*nhh{1},t2pmod2{1},Attmod2{1},levs(:));
Eout=Eout1;rs=rs1;js=js1;

levs=cat(1,clevs(Nv(1)+1:Nv(1)+Nv(2)),-clevs(Nv(1)+1:Nv(1)+Nv(2)));
[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{2}-.002*nhh{2},t2pmod2{2},Attmod2{2},levs(:));
Eout=Eout+Eout1;rs=cat(1,rs,rs1);js=cat(1,js,js1);

[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{2}+.002*nhh{2},t2pmod2{2},Attmod2{2},levs(:));
Eout=Eout+Eout1;rs=cat(1,rs,rs1);js=cat(1,js,js1);

levs=cat(1,clevs(Nv(1)+Nv(2)+1:end),-clevs(Nv(1)+Nv(2)+1:end));
[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{3}-.002*nhh{3},t2pmod2{3},Attmod2{3},levs(:));
Eout=Eout+Eout1;rs=cat(1,rs,rs1);js=cat(1,js,js1);

Eout=Eout*Jval;
c=norm(Eideal(:)-Eout(:))/norm(Eideal(:));%+real(computeinductance(rs,js))*10^5/5;
c
end