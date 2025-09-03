
function [p,t2p]=addcube(p,t2p,locz)
np2=numel(p)/3;
p(end+1:end+8,1)=.01*[0,1,0,1,0,1,0,1];
p(np2+1:np2+8,2)=.01*[0,0,1,1,0,0,1,1];
p(np2+1:np2+8,3)=.01*[0,0,0,0,1,1,1,1]+locz;
t2p(end+1,:)=[np2+1;np2+2;np2+3];
t2p(end+1,:)=[np2+2;np2+4;np2+3];
t2p(end+1,:)=[np2+7;np2+6;np2+5];
t2p(end+1,:)=[np2+6;np2+7;np2+8];

t2p(end+1,:)=[np2+1;np2+2;np2+5];
t2p(end+1,:)=[np2+5;np2+6;np2+2];
t2p(end+1,:)=[np2+3;np2+4;np2+7];
t2p(end+1,:)=[np2+7;np2+8;np2+4];

t2p(end+1,:)=[np2+1;np2+3;np2+5];
t2p(end+1,:)=[np2+5;np2+7;np2+3];
t2p(end+1,:)=[np2+2;np2+4;np2+6];
t2p(end+1,:)=[np2+6;np2+8;np2+4];
