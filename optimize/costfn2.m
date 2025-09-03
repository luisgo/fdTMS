function c=costfn2(levs,Eout,psym,t2preg,Att1,nx,ny,nz)
levs=cat(1,levs(:),-levs(:));
[rs,js,Eout2,~,~,~,~,~,~,~,boo]=coilfieldsp(psym,t2preg,Att1,levs(:),nx,ny,nz);
c=norm(Eout(:)-Eout2(:))/norm(Eout(:))+real(computeinductance(rs,js))*10^5/5;

end


