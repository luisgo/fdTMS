function c=costfnall(clevs,Eideal,pmod2,t2pmod2,Attmod2,Jval,Nv,nhh)

levs=cat(1,clevs(1:Nv(1)),-clevs(1:Nv(1)));
[rs1,js1,Eout,p,t2p,clen,...
    xcoilbig,ycoilbig,zcoilbig,zzcoilbig]=coilfields(pmod2{1},t2pmod2{1},Attmod2{1},levs(:));


levs=cat(1,clevs(Nv(1)+1:Nv(1)+Nv(2)),-clevs(Nv(1)+1:Nv(1)+Nv(2)));
[rs2,js2,Eout,p,t2p,clen,...
    xcoilfig8,ycoilfig8,zcoilfig8,zzcoilfig8]=coilfields(pmod2{2},t2pmod2{2},Attmod2{2},levs(:));

levs=cat(1,clevs(Nv(1)+Nv(2)+1:end),-clevs(Nv(1)+Nv(2)+1:end));
[rs3,js3,~,p,t2p,clen,...
    xcoilside,ycoilside,zcoilside,zzcoilside]=coilfields(pmod2{3},t2pmod2{3},Attmod2{3},levs(:));

[Eout,rs,js]=optimseparate(xcoilbig,ycoilbig,zcoilbig,...
          xcoilfig8,ycoilfig8,zcoilfig8,...
          xcoilside,ycoilside,zcoilside);
Eout=Eout*Jval;
c=norm(Eideal(:)-Eout(:))/norm(Eideal(:));%+real(computeinductance(rs,js))*10^5/5;
c
end