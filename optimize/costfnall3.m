function c=costfnall3(clevs,Eideal,pmod2,t2pmod2,Attmod2,Jval,Nv,nhh)

levs=cat(1,clevs(1:Nv(1)),-clevs(1:Nv(1)));
[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{1}-.0025*nhh{1},t2pmod2{1},Attmod2{1},levs(:));
Eout=Eout1;rs=rs1;js=js1;

levs=cat(1,clevs(Nv(1)+1:Nv(1)+Nv(2)),-clevs(Nv(1)+1:Nv(1)+Nv(2)));
[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{2}-.0025*nhh{2},t2pmod2{2},Attmod2{2},levs(:));
Eout=Eout+Eout1;rs=cat(1,rs,rs1);js=cat(1,js,js1);

[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{2}+.0005*nhh{2},t2pmod2{2},Attmod2{2},levs(:));
Eout=Eout+Eout1;rs=cat(1,rs,rs1);js=cat(1,js,js1);

[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{2}+.0035*nhh{2},t2pmod2{2},Attmod2{2},levs(:));
Eout=Eout+Eout1;rs=cat(1,rs,rs1);js=cat(1,js,js1);

levs=cat(1,clevs(Nv(1)+Nv(2)+1:end),-clevs(Nv(1)+Nv(2)+1:end));
[rs1,js1,Eout1,~,~,~,~,~,~,~]=coilfields(pmod2{3}-.0025*nhh{3},t2pmod2{3},Attmod2{3},levs(:));
Eout=Eout+Eout1;rs=cat(1,rs,rs1);js=cat(1,js,js1);
Eout=Eout*Jval;
Lest=real(computeinductance(rs,js))
West=Lest*Jval^2/2
c=norm(Eideal(:)-Eout(:))/norm(Eideal(:))+.001*West/200;

c
end