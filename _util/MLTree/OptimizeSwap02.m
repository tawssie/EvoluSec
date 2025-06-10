function SWAP=OptimizeSwap02(TreeLikelihoodObjFun,Model,CMatrix,ShortestPath,NodeIndex,LSet,LSetLeft,LSetRight,P)
% This function is to optimize single edge by all pre-calculated information from GetTreeLikelihood03.
% WITH alternative 2 topologies from SWAPS
%
% Input:
%     TreeLikelihoodObjFun
%                   <= obj function to optimize single edge by all pre-calculated information as follows
%     Model         <= evolutionary model
%     CMatrix       <= known tree topology
%     ShortestPath  <= pre-calculated shortest path from "NodeIndex" to the root
%     LSet, LSetLeft, LSetRight <= from previous CondL structure, with conditional likelihood
%     P             <= square matrix of cells with all probability matrix
% Output:
%     SWAPList      <= a cell list, each cell contains two optimization information for 2 swaps
% 
% It's NOT neccessary to update full condiontal likelihood (LSet*)
% cause only Cond Ls along the shortest path from node to root needs to re-calculated by optimize edge
% the point is to guarantee information aroung optimized edge as follows:
% 1. alternative topology AlterM1, AlterM2 are correct
% 2. maintain relationships between Current Node and New brother from their parent
%     (cause parent's left child index < right child index, but it will be violated after swaps)
% 3. to resolve 2., reassign parent's LSetLeft (or LSetRight) by Current Node and New brother
% 4. LSet of Current Node (NodeIndex) is UPDATED according to alternative topology (LSet=child.*old_brother)



ChildList=GetChild(CMatrix,NodeIndex);
% Smaller index (first element) is always the left child
ChildLeftIndex =ChildList(1);
ChildRightIndex=ChildList(2); 

ParentIndex =GetParent (CMatrix,NodeIndex);
BrotherIndex=GetBrother(CMatrix,NodeIndex); %%% brother is the "neighbor" we're going to use


% Alternative topology (swap) 1, exchange ChildLeft and Brother
%     i.  update LSet for Current Node
%     ii. update Parent Node's LSetLeft (or LSetRight)
%     iii.optimize edge length for Current Node

% Get a copy of all conditional likelihood, modify them for later optimization
AlterLSet1     =LSet;
AlterLSetLeft1 =LSetLeft;
AlterLSetRight1=LSetRight;

% Do NNI, AlterCM1 is the new topology, exchanging ChildLeft and Brother
AlterCM1=CMatrix;
AlterCM1(ParentIndex,BrotherIndex  )=0;
AlterCM1(NodeIndex  ,BrotherIndex  )=1;
AlterCM1(NodeIndex  ,ChildLeftIndex)=0;
AlterCM1(ParentIndex,ChildLeftIndex)=1;

% Get a copy of all transition matrix, modify them according to new topology for later optimization
AlterP1=P;
AlterP1(NodeIndex  ,BrotherIndex  )=P(ParentIndex,BrotherIndex  );
AlterP1(ParentIndex,ChildLeftIndex)=P(NodeIndex  ,ChildLeftIndex);
AlterP1(ParentIndex,BrotherIndex  )={[]};
AlterP1(NodeIndex  ,ChildLeftIndex)={[]};


% 1. Update LSet of NodeIndex:
%   1.1 According to old topology, test whether brother is left child of parent.
%    Retrive brother's cond L (TempL) from its parent
if IsLeftChild(CMatrix,ParentIndex,BrotherIndex)
    TempL=LSetLeft (ParentIndex,:);
else
    TempL=LSetRight(ParentIndex,:);
end

%   1.2. re-calculate LSet for NodeIndex
for i=1:size(AlterLSet1,2)
    AlterLSet1{NodeIndex,i}=TempL{i}.*LSetRight{NodeIndex,i};
end


% 2. Update Parent Node's LSetLeft (or LSetRight)
%   1.1 Since AlterCM1 is to exchange ChildLeft with Brother
%       retrieve ChildLeft's cond L from Current Node
TempL=LSetLeft(NodeIndex,:);

%   1.2 According to New topology, test whether Current Node is left child of parent.
%       if current is left child, brother is right child, so TempL should be assigned to right  
if IsLeftChild(AlterCM1,ParentIndex,NodeIndex)
    % in new topology, NodeIndex is left child of parent
    % so update right condL of parent
    % cause edge with NodeIndex will be optimized
    % but right child of parent remain the same
    AlterLSetRight1(ParentIndex,:)=TempL;% warning('case 1');
else
    AlterLSetLeft1 (ParentIndex,:)=TempL; 
end


% optimize edge length for current node
[D1,TreeLogL1,~,~,~]=BrentOptAlgo4Tree(TreeLikelihoodObjFun...
    ,0,6,Model,AlterCM1,ShortestPath,AlterLSet1,AlterLSetLeft1,AlterLSetRight1,AlterP1);

%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%


% alternative topology (swap) 2
% Get a copy of all conditional likelihood, modify them for later optimization
% Do NNI, AlterCM2 is the new topology, exchanging ChildLeft and Brother
% Get a copy of all transition matrix, modify them according to new topology for later optimization
AlterLSet2=LSet;
AlterLSetLeft2=LSetLeft;
AlterLSetRight2=LSetRight;

AlterCM2=CMatrix;
AlterCM2(ParentIndex,BrotherIndex   )=0;
AlterCM2(NodeIndex  ,BrotherIndex   )=1;
AlterCM2(NodeIndex  ,ChildRightIndex)=0;
AlterCM2(ParentIndex,ChildRightIndex)=1;

AlterP2=P;
AlterP2(NodeIndex  ,BrotherIndex   )=P(ParentIndex,BrotherIndex   );
AlterP2(ParentIndex,ChildRightIndex)=P(NodeIndex  ,ChildRightIndex);
AlterP2(ParentIndex,BrotherIndex   )={[]};
AlterP2(NodeIndex  ,ChildRightIndex)={[]};


if IsLeftChild(CMatrix,ParentIndex,BrotherIndex)
    TempL=LSetLeft (ParentIndex,:);
else
    TempL=LSetRight(ParentIndex,:);
end

%   re-calculate LSet for NodeIndex
for i=1:size(AlterLSet2,2)
    AlterLSet2{NodeIndex,i}=LSetLeft{NodeIndex,i}.*TempL{i};
end

TempL=LSetRight(NodeIndex,:);

if IsLeftChild(AlterCM2,ParentIndex,NodeIndex)
    AlterLSetRight2(ParentIndex,:)=TempL;    
else
    AlterLSetLeft2 (ParentIndex,:)=TempL;
end


% optimize edge length for current node
[D2,TreeLogL2,~,~,~]=BrentOptAlgo4Tree(TreeLikelihoodObjFun...
    ,0,6,Model,AlterCM2,ShortestPath,AlterLSet2,AlterLSetLeft2,AlterLSetRight2,AlterP2);


% SWAP is a record with optimized edge information corresponding to 2 possible swaps
SWAP=...
[D1 TreeLogL1 NodeIndex ParentIndex ChildLeftIndex  BrotherIndex
 D2 TreeLogL2 NodeIndex ParentIndex ChildRightIndex BrotherIndex];

    
    
function ParentIndex=GetParent(M,NodeIndex)
% ParentIndex=[], if no parent => the root
ParentIndex=find(M(:,NodeIndex));


function BrotherIndex=GetBrother(M,NodeIndex)
% BrotherIndex=[],if no brother => the root

% this is children of parent (1 is NodeIndex, another is brother)
Child=find(M(M(:,NodeIndex),:)); 


% use numel to test both root and single leaf condition
% numel([]) == 0, so numel can both test isempty and zero length
% numel(Child) == 0 , NodeIndex == the root

if numel(Child)==0
    % it happens when NodeIndex is the root
    BrotherIndex=[];
else
    if Child(1)==NodeIndex
        BrotherIndex=Child(2);
    else
        BrotherIndex=Child(1);
    end

end

function ChildIndexList=GetChild(M,NodeIndex)

ChildIndexList=find(M(NodeIndex,:));

function i=IsLeftChild(M,NodeIndex,ChildIndex)

ChildList=GetChild(M,NodeIndex);

i=false;
if ChildIndex==ChildList(1)
    i=true;
end
