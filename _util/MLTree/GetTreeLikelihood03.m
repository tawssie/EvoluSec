function [TreeLogL,CondL,P]=GetTreeLikelihood03(TipLikelihood,MSALabel,Model,CMatrix,TraverseOrder,D)
% This is the function to "calculate" tree **LOG** likelihood by
%     1. Topology <= from CMatrix
%     2. Pre-defined likelihood for each tip <= TipLikelihood (this is NOT in log scale)
%     3. Evolutionary model e.g., JTT, Dayhoff, SS7Label <= Model
%     4. Pre-calculated "reverse" traverse route as DFS from root <= TraverseOrder (save time in purpose)
%     5. Evolutionary distance <= D, D(i) is the evolutionary distance from node_i's parent to node_i
%
%
% The function returns 3 object
%     1, TreeLogL => loglikelihood of a whole tree
%     2. CondL    => all conditional likelihood. a struct with LSet, LSetRight, LSetLeft
%            LSet      (conditional likelihood at node i)
%            LSetLeft  (conditional likelihood at node i's left  child != left  child's cond likelihood)
%            LSetRight (conditional likelihood at node i's right child != right child's cond likelihood)
%            LSet{5,1} == (LSetRight{5,1}.*LSetRight{5,1}) node 5 at column 1 of MSA
%     3. P        => square cell matrix same dimension as CMatrix,
%            P{i,j} is the prob transition matrix from node i to node j
%
% TipLikelihood differs in different evolutionary model, 
%     details in GetTipsLikelihood4AminoAcid and GetTipsLikelihood4SS


% Calculates transition probability matrix for all edges (nodes) as P (cell matrix)
% GetExpMByTime return as a cell list when numel(t) > 1
P=cell(size(CMatrix));
P(CMatrix==1)=GetExpMByTime(Model.V,Model.D,D);

% LSet records all conditional likelihood = numel(total tree nodes) X numel(MSA columns)
%   "row ordering" followed by node number
%   NOT traversal order, i.e., row 5 = node 5's conditional L
LSet     =cell(size(CMatrix,1),size(MSALabel,2));
LSetLeft =cell(size(CMatrix,1),size(MSALabel,2));
LSetRight=cell(size(CMatrix,1),size(MSALabel,2));


% Pre calculate child list for all nodes (internal + tips) to speed up
ChildList=cell(numel(TraverseOrder),1);
for i=1:numel(TraverseOrder)
    %since i use "find" to define child list
    %ChildList{i} is 1*2 row vector
    %and Left Child is defined as smaller index (fist element of list)
    ChildList{i}=find(CMatrix(TraverseOrder(i),:));
end


% 1. This version introduces GapL to avoid precision loss.
% 2. Also rearrange i, c position to improve speed
GapL=ones(numel(Model.F),1); % Gap vectors's dimension == dimension of Model.F
for i=1:numel(TraverseOrder)
    if isempty(ChildList{i})
        for c=1:size(MSALabel,2)
            LSet{TraverseOrder(i),c}=TipLikelihood(:,MSALabel(TraverseOrder(i),c));
        end
    else
        PLeftChild=P{TraverseOrder(i),ChildList{i}(1)};
        PRightChild=P{TraverseOrder(i),ChildList{i}(2)};
        for c=1:size(MSALabel,2)
            if isequal(LSet{ChildList{i}(1),c},GapL)
                LSetLeft {TraverseOrder(i),c}=GapL;
            else
                LSetLeft {TraverseOrder(i),c}=PLeftChild *LSet{ChildList{i}(1),c};
            end
            
            if isequal(LSet{ChildList{i}(2),c},GapL)
                LSetRight{TraverseOrder(i),c}=GapL;
            else
                LSetRight{TraverseOrder(i),c}=PRightChild*LSet{ChildList{i}(2),c};
            end

            LSet{TraverseOrder(i),c}=LSetLeft {TraverseOrder(i),c}.* LSetRight{TraverseOrder(i),c};
            
        end
    end
end


% Below is matrix style operation to retrieve tree log likelihood
F=Model.F(:)'; % F is 1*20 matrix of equlibrium frequency
               % Project all TREE likelihood vector to F in the next step

% The end of TraverseOrder is the root so below
% is a matrix contains all TREE likelihood vector
% [LSet{TraverseOrder(end),:}]; 
TreeLogL=sum( log(F*[LSet{TraverseOrder(end),:}]) );

CondL.LSet=LSet;
CondL.LSetLeft=LSetLeft;
CondL.LSetRight=LSetRight;

