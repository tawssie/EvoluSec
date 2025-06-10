function [AlterM,M,ID,Dist]=NNInterchange03(PhyTreeObject,ID,Dist)
% [AlterM,M,ID,Dist]=NNInterchange03(tree05);
% [AlterM]=NNInterchange03(M,ID,Dist);
%
% AlterM is a cell matrix of NumOfNode bt 2
% AlterM(i,1), AlterM(i,2) are two alternative topology
% at the edge connecting to node i
% if node i is root or child, Alter(i,1), Alter(i,2) will be empty {[]}
%
% nearest neighbor interchanges
% the edge is specified by NodeIndex
% , and edge attaches to parent of NodeIndex
%
% generate the other two alternative topology according to NodeInde
% each alternative has its tree likelihood
%
% a bifurcated phylogenetic tree, each internal nodes MUST have two children



if nargin==1&&isa(PhyTreeObject,'phytree')
    [M,ID,Dist]=getmatrix(PhyTreeObject);
    M=full(M);
elseif nargin==1&&isa(PhyTreeObject,'logical')
    M=PhyTreeObject;
    M=full(M);
    ID=[];
    Dist=[];
elseif nargin==3
    M=PhyTreeObject;
else
    error('wrong input argument');
end

AlterM=cell(size(M,1),2);


for i=find(any(M,2)&any(M,1)')'
%%% Go through internal nodes

    MM=GetNNIByNodeIndex(M,i);
    
    if ~isempty(MM)
        AlterM(i,:)=MM;
    end
end


function MM=GetNNIByNodeIndex(M,NodeIndex)
% 1. Only works for bifurcated phylogenetic tree
% 2. MM = [], if NodeIndex is root


MM=[];

ChildList=find(M(NodeIndex,:));

NumOfChild=numel(ChildList);

if NumOfChild==2
    
    MM=cell(2,1);

    ParentIndex =GetParent (M,NodeIndex);
    BrotherIndex=GetBrother(M,NodeIndex);
    

    MM1=M;
    MM1(ParentIndex,BrotherIndex)=0;
    MM1(NodeIndex  ,BrotherIndex)=1;

    MM1(NodeIndex  ,ChildList(1))=0;
    MM1(ParentIndex,ChildList(1))=1;

    MM{1}=MM1;
    
    MM2=M;
    MM2(ParentIndex,BrotherIndex)=0;
    MM2(NodeIndex  ,BrotherIndex)=1;

    MM2(NodeIndex  ,ChildList(2))=0;
    MM2(ParentIndex,ChildList(2))=1;

    MM{2}=MM2;

end



function ParentIndex=GetParent(M,NodeIndex)
% ParentIndex=[], if no parent => the root
ParentIndex=find(M(:,NodeIndex));


function BrotherIndex=GetBrother(M,NodeIndex)
% BrotherIndex=[],if no brother => the root

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


