function [xcoil,ycoil,zcoil]=resamplecoil(xcoilbig,ycoilbig,zcoilbig,ids);
ct=1;
for i=ids
xcoil{ct}=xcoilbig{i};
ycoil{ct}=ycoilbig{i};
zcoil{ct}=zcoilbig{i};
ct=ct+1;
      end
clen=0;
resolu=.0001;
for i=1:numel(xcoil)
    len=sqrt((xcoil{i}(2:1:end)-xcoil{i}(1:end-1)).^2+...
             (ycoil{i}(2:1:end)-ycoil{i}(1:end-1)).^2+...
             (zcoil{i}(2:1:end)-zcoil{i}(1:end-1)).^2);    
         
     len=cumsum([0;len(:)]);
     clen=clen+len(end);
     N=ceil(len(end)/resolu);
     locs=0:len(end)/N:len(end);
     xcoil{i}=interp1(len,xcoil{i},locs,'linear');
     ycoil{i}=interp1(len,ycoil{i},locs,'linear');
     zcoil{i}=interp1(len,zcoil{i},locs,'linear');    
end