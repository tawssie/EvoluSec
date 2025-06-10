function [NewCMatrix,NewD,SWAPIndex]=ImproveTree04(TreeLogL,CondL,P,Model,CMatrix,D,DoParfor)
% Output:
%     NewCMatrix = updated topology
%     NewD       = updated distance list (this will be further modified by Lambda)
%     SWAPIndex  = used to determine which kind of lambda will apply on NewD
% Input:
%     TreeL   = current tree loglikelihood
%     CondL   = struct with all conditional likelihood, LSet, LSetRight, LSetLeft
%     P       = cell matrix, each cell contains a probability matrix corresponding to edge length
%               P{i,j}== probability matrix from parent i to child j with distance between i and j
%     Model   = evolutionary model, could be amino acid or secondary structure
%     CMatrix = connected matrix
%     TraverseOrder = reverse DFS from "node" that "node's edge" is going to be optimized.
%                     first element = tips; last element = pseudo root
% 
%     D       = distance list, D(i) = edge length at node_i
% 
%     since
%         P{i,j}== probability matrix from parent i to child j
%     and
%         D(j)== edge attached on the node j
%     therefore,
%         GetExpMByTime(Model.V,Model.D,D(j)) == P{i,j}


% Retrieve all conditional probabilities
LSet     =CondL.LSet;
LSetLeft =CondL.LSetLeft;
LSetRight=CondL.LSetRight;

% Retrieve all shortest path as a list for all nodes except for root
ShortestPathList=GetShortestPathListAllNode2Root(CMatrix);


% 1.  Optimize edge length without SWAPPING
%     Individually optimize edge length (evolutionary distance from parent of node_i to current node_i)
%     ShortestPathList{i} is the traverse route from node_i to the root
%         , i.e., conditional likelihood must be updated for all edges within the path
%     CAUTION,
%         TreeLogL is the "additive inverse, x=-x" because we obtain maximum value by "searching the minimum"

D_noSWAP=zeros(numel(ShortestPathList),1);
TreeLogL_noSWAP=zeros(numel(ShortestPathList),1);

if DoParfor
    parfor i=1:numel(ShortestPathList)  % parallel computing
        [D_noSWAP(i),TreeLogL_noSWAP(i),~,~,~]=BrentOptAlgo4Tree(@TreeLikelihoodObjFun...
            ,0,6,Model,CMatrix,ShortestPathList{i},LSet,LSetLeft,LSetRight,P);
    end
else
    for i=1:numel(ShortestPathList)
        [D_noSWAP(i),TreeLogL_noSWAP(i),~,~,~]=BrentOptAlgo4Tree(@TreeLikelihoodObjFun...
            ,0,6,Model,CMatrix,ShortestPathList{i},LSet,LSetLeft,LSetRight,P);
    end
end

% Test whether TreeLogL is a complex or a NaN.
% It shouldn't happen but I keep it for additional checking.
if any(~isreal(TreeLogL_noSWAP)), warning('Optimized Tree log Likelihood by single edge is not a real number'); end
if any(  isnan(TreeLogL_noSWAP)), warning('Optimized Tree log Likelihood by single edge is NaN'); end


% 2.  Optimize edge length with 2 possible nearest-neighbor interexchange
%     Individually optimize edge length with Nearest Neighbor Interexchange
%         (evolutionary distance from parent of node_i to current node_i)
%     SWAPList == cell list to store all potential swappings with their
%         1. Optimized edge length
%         2. TreeLogLikelihood (additive inverse)
%         3. NodeIndex
%         4. ParentIndex (parent of NodeIndex)
%         5. ChildLeft or ChildRight (two possible swappings for NodeIndex)
%         6. BrotherIndex (brother/nearest neighbor of NodeIndex)


SWAPList=cell(numel(ShortestPathList),1);


if DoParfor
    parfor NodeIndex=find(any(CMatrix,2)'&any(CMatrix,1)) % parallel computing
        SWAPList{NodeIndex}=OptimizeSwap02...
            (@TreeLikelihoodObjFun...
            ,Model,CMatrix,ShortestPathList{NodeIndex},NodeIndex,LSet,LSetLeft,LSetRight,P);
    end
else
    for NodeIndex=find(any(CMatrix,2)'&any(CMatrix,1))
        SWAPList{NodeIndex}=OptimizeSwap02...
            (@TreeLikelihoodObjFun...
            ,Model,CMatrix,ShortestPathList{NodeIndex},NodeIndex,LSet,LSetLeft,LSetRight,P);
    end
end


% 3.  Combine results from the above two optimizations
%     {D, TreeLogL, CMatrix},they are from previous stage
%         D        = edge length for all nodes except root in the tree
%         TreeLogL = log tree likelihood in additive inverse
%         CMatrix  = tree topology
%     D_noSWAP, TreeLogL_noSWAP <= optimization without SWAP
%     SWAPList <= cell list containing all information from swapping


[NewCMatrix,NewD,SWAPIndex]=CombineOptimization02(D,TreeLogL,CMatrix,D_noSWAP,TreeLogL_noSWAP,SWAPList);






function ExtraResult=TreeLikelihoodObjFun(D,Model,CMatrix,TraverseOrder,LSet,LSetLeft,LSetRight,P)
% This function is similar to GetTreeLikelihood03 (calculate tree loglikelihood with MSA, topology, distance)
% BUT here is used to optimize single edge by all pre-calculated information from GetTreeLikelihood03
% ***CAUTION*** D is moved to the first parameter in this function
%     in order to mathch ExtraResult parameter within fminbnd_JSLai03
%
% Details in GetTreeLikelihood03.m
% 
% D             <= unknown single edge length (distance), to be optimized by fminbnd_JSLai03
% Model         <= evolutionary model
% CMatrix       <= known tree topology
% TraverseOrder <= pre-calculated shortest path from "node attached to D" to the root
% LSet, LSetLeft, LSetRight <= from previous CondL structure, with conditional likelihood
% P             <= square matrix of cells with all probability matrix


NewP=GetExpMByTime(Model.V,Model.D,D);

% Update transition probability matrix directly
P{CMatrix(:,TraverseOrder(1)),TraverseOrder(1)}=NewP;


% Pre calculate child list for all nodes (internal + tips) to speed up
ChildList=cell(numel(TraverseOrder),1);
for i=1:numel(TraverseOrder)
    %Since i use "find" to define child list
    %ChildList{i} is 1*2 row vector
    %And "Left Child" is defined as smaller index (fist element of list)
    ChildList{i}=find(CMatrix(TraverseOrder(i),:));
end


% 1. No parallel computing here, parallel computing in "searching topology"
% 2. This version reorders left child testing, in order to avoid computation overhead
% 3. Also use GapL to assign probability directly for Gap node, to avoid precision loss
%    P_any_time * GapL == GapL because it just sum by columns of P_any_time, should be "one vector"
GapL=ones(numel(Model.F),1); % Gap vectors's dimension == dimension of Model.F 


% No need to update the 1st element of traverse order, its conditional likelihood is fine
% so iterator "i" starts from  2
for i=2:numel(TraverseOrder)    
    % Current node is TraverseOrder(2)
    % TraverseOrder(i-1) is child of TraverseOrder(i) because
    % TraverOrder is the shortest path from node to the root
    
    PLeftChild=P{TraverseOrder(i),ChildList{i}(1)};
    PRightChild=P{TraverseOrder(i),ChildList{i}(2)};
    
    % Test Child == left child    
    if TraverseOrder(i-1)==ChildList{i}(1)
        for c=1:size(LSet,2) % go throught all columns of MSA
            % Update Conditional likelihood:
            %       Update 1. LSetLeft/LSetRight (depends on child)
 
            % special treatment for gap character
            if isequal(LSet{ChildList{i}(1),c},GapL)
                % assign probabilities "1" to gap observed likelihood to avoid precision loss
                LSetLeft {TraverseOrder(i),c}=GapL;
            else                
                LSetLeft {TraverseOrder(i),c}=PLeftChild *LSet{ChildList{i}(1),c};
            end
            %       Update 2. LSet (node), 
            LSet{TraverseOrder(i),c}=LSetLeft{TraverseOrder(i),c}.* LSetRight{TraverseOrder(i),c};
        end
    else % Otherwise, update right child ChildList{i}(2)
        for c=1:size(LSet,2)
            if isequal(LSet{ChildList{i}(2),c},GapL)
                LSetRight{TraverseOrder(i),c}=GapL;
            else
                LSetRight{TraverseOrder(i),c}=PRightChild*LSet{ChildList{i}(2),c};
            end
            LSet{TraverseOrder(i),c}=LSetLeft{TraverseOrder(i),c}.* LSetRight{TraverseOrder(i),c};
        end
    end
end


% Below is matrix style operation to retrieve tree log likelihood
F=Model.F(:)'; % F is 1*20 matrix of equlibrium frequency
               % Project all TREE likelihood vector to F in the next step
TreeLogL=sum( log(F*[LSet{TraverseOrder(end),:}]) );

% Additive inverse, in order to apply "minimum" searching algorithm
TreeLogL=-TreeLogL;

ExtraResult={TreeLogL LSet LSetLeft LSetRight P};



function ShortestPathList=GetShortestPathListAllNode2Root(CMatrix)
% This is a funcion to get shortest path from node to the root
% ShortestPathList{i} == 1*n, row vector represents traverse order from
%                        node_i to the root (RootIndex), 

RootIndex=find(~any(CMatrix,1));

M=sparse(CMatrix|CMatrix');

NumOfNode=size(CMatrix,1)-1; % # of internal + # of tips, EXCEPT ROOT
ShortestPathList=cell(NumOfNode,1);

for i=1:NumOfNode
    % Record all path from root to all nodes (both tips and internal)
    [~,ShortestPathList{i}]=graphshortestpath(M,RootIndex,i);
    ShortestPathList{i}=fliplr(ShortestPathList{i});
end

