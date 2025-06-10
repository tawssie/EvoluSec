function tt=CreatePhytreeFromBGraph(CMatrix,ID,Dist)
CMatrix=sparse(CMatrix);

TotalLeaves=full(sum(~any(CMatrix,2)));
RootIndex=find(~any(CMatrix,1));

z=flipud(graphtraverse(CMatrix,RootIndex,'method','BFS')');


zz=z;

TotalNodes=size(CMatrix,1);

zz(zz<=TotalLeaves)=1:TotalLeaves;
zz(zz>TotalLeaves)=TotalLeaves+1:TotalNodes;

c=[zz z];
c=sortrows(c,1);

ID2=ID(c(:,2));
Dist2=Dist(c(:,2));

k=reshape(zz(1:end-1),2,TotalLeaves-1);
tt=phytree(k',Dist2,ID2);

