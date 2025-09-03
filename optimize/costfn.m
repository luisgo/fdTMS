function c=costfn(levs,Eout,psym,t2preg,Att1)
levs=cat(1,levs(:),-levs(:));
[rs,js,Eout1,~,~,~,~,~,~,~,boo]=coilfields(psym{1},t2preg{1},Att1{1},levs(:));
c=norm(Eout(:)-Eout2(:))/norm(Eout(:))+real(computeinductance(rs,js))*10^5/5;

end