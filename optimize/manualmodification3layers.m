function wireth=manualmodification3layers(wireth1,wireth2,wid1,wid2,rs);
wireth=wireth1*ones(size(rs(:,1)));
ct=0;
for j=[1 2 3 4 5 6 7 8 10 11 21]
for i=wid1(j):wid2(j)
if ct<100
wireth(i)=wireth1;
if ct==99 && i>=wid2(j)-120
    break
end
elseif ct<200
wireth(i)=wireth2;
end
ct=ct+1;
if ct==200
ct=0;
end
end
ct=0;
end
