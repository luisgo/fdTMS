
function getmodes(i)
if i==1
p=load('hsphere.nodes');
p=p(:,end-2:end);
t2p=load('hsphere.triag');
t2p=t2p(:,end-2:end);
trisurf(t2p,p(:,1),p(:,2),p(:,3))
addpath('../generateharmonics/');
generateharmonicsinputmeshtri('filamen_halfsphere',1,p,t2p,0,30);
elseif i==2
p=load('hsphere.nodes');
p=p(:,end-2:end);
t2p=load('hsphere.triag');
t2p=t2p(:,end-2:end);
trisurf(t2p,p(:,1),p(:,2),p(:,3))
addpath('../generateharmonics/');
generateharmonicsinputmeshtri('thick_halfsphere',1,.93*p/.9,t2p,0,30);

elseif i==3
p=load('spherem.nodes');
p=p(:,end-2:end);
t2p=load('spherem.triag');
t2p=t2p(:,end-2:end);
trisurf(t2p,p(:,1),p(:,2),p(:,3))
generateharmonicsinputmeshtri('filamen_sphere',1,p,t2p,0,30);

elseif i==4
p=load('spherem.nodes');
p=p(:,end-2:end);
t2p=load('spherem.triag');
t2p=t2p(:,end-2:end);
trisurf(t2p,p(:,1),p(:,2),p(:,3))
generateharmonicsinputmeshtri('thick_sphere',1,.93*p/.9,t2p,0,30);

elseif i==5
N=70-1;
[X,Y]=ndgrid(-.16:.16/N:0,0:.16/N:.16);
t2p=delaunay(X(:),Y(:))
p(:,1)=X(:);
p(:,2)=Y(:);
p(:,3)=.09;

generateharmonicsinputmeshtri('filamen_square',1,p,t2p,0,30);

elseif i==6
N=70-1;
[X,Y]=ndgrid(-.16:.16/N:0,0:.16/N:.16);
t2p=delaunay(X(:),Y(:))
p(:,1)=X(:);
p(:,2)=Y(:);
p(:,3)=.093;
generateharmonicsinputmeshtri('thick_square',1,.93*p/.9,t2p,0,30);

end
