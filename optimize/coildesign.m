[xcoil,ycoil,zcoil,ff,Fsym1,Fsym2,Att1,psym,t2psym...
    ]=getcoilwindings(harmfile,designfile,1,2,10,16);
close all
for i=1:numel(xcoil)
   plot(xcoil{i},ycoil{i})
   hold on
   text(xcoil{i}(1),ycoil{i}(1),num2str(i)) 
end