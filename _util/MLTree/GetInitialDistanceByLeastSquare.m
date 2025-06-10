function D_init=GetInitialDistanceByLeastSquare(CMatrix,PairwiseDist)
% Use "lsqnonneg" function to solve least-squares constraints with edge length >=0 
%
% Folling MOLPHY's definition
% MOLPHY version 2.3: programs for molecular phylogenetics based on maximum likelihood
%     A = topology matrix A of MOLPHY
%     PairwiseDist = D of MOLPHY
%     D_init = estimated t of MOLPHY, node 1~n distance
%     D_init(5) = branch length connects to node 5

M=sparse(CMatrix|CMatrix'); % make it symmetric and sparse

NumOfTips=sum(any(CMatrix,2)==0); % count # of leaves

% Determine shortest path between every pair of tip nodes with the tree.
[y,x]=find(tril(true(NumOfTips),-1));   
GraphPath=cell(numel(x),1);

for i=1:numel(GraphPath)
    [~,GraphPath{i}]=graphshortestpath(M,x(i),y(i));
end


% Determine topology matrix
A=zeros(nchoosek(NumOfTips,2),NumOfTips+(NumOfTips-1));


for i=1:numel(GraphPath)
    P=GraphPath{i}; % P is the path, from x(i) => y(i)
                    % the path records NODES of tree
                    % so use min(.) to determin edge from two ADJACENT nodes
    for j=1:numel(P)-1
        A(i,min(P(j:j+1)))=1;
    end
end

D_init=lsqnonneg(A,PairwiseDist);%%% distances are constrained >=0 



