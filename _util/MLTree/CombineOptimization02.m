function [NewCMatrix,NewD,SWAPIndex]=CombineOptimization02(D,TreeLogL,CMatrix,D_noSWAP,TreeLogL_noSWAP,SWAPList)
% This function combines optimizations of all single edges from non-swapping and swapping results
% And returns new topology, optimized edge length, and index to represent swapping.
% Following PHYML's description to make sure
% 1. swapping's priority > non-swapping
% 2. only keep distances that can improve tree log likelihood
%
% Input:
%     at previous stage
%     D        = edge length, to decide whether optimize distance should be updated
%     TreeLogL = tree log likelihood
%     CMatrix  = tree topology
%     
%     at current stage
%     D_noSWAP = optimized edge length for all nodes except root
%     TreeLogL_noSWAP
%              = tree log likelihood corresponding to D_noSWAP
%     SWAPList = cell list contains information from NNI swapping
% Output:
%     NewCMatrix = new topology (swapping leads to new tree)
%     NewD       = "potential" optimized edge length, will be further scaled by lambda
%     SWAPIndex  = index of {0,1,2} to determine lambda for each edge in the later




% Make two tables from both swap and non-swap results
% I combine them all together for further comparsion with tree log likelihood
% Table_noSWAP Table_SWAP share same format
%     1. edge_length
%     2. tree_logL
%     3. current node index
%     4. parent node index
%     5. child left or child right
%     6. brother node index
% if no swap, current node == parent node == child node ...
%
% Besides, make sure TreeLogL_noSWAP takes additive inverse
% cause i use "searching minimum" to maximize tree likelihood.

Table_noSWAP=[D_noSWAP -TreeLogL_noSWAP repmat([1:numel(D)]',1,4)];
Table_SWAP=vertcat(SWAPList{:});
Table_SWAP(:,2)=-Table_SWAP(:,2);

T=CombineTable(Table_noSWAP,Table_SWAP,TreeLogL);


% indicator if current node index != parent node index
IsSWAP=T(:,3)~=T(:,4); 
SWAPRecords   =T( IsSWAP,:);
NonSWAPRecords=T(~IsSWAP,:);

NewD=D; %important to keep old D
% update distance (:,1) by current node index (:,3)
NewD(NonSWAPRecords(:,3))=NonSWAPRecords(:,1);
NewD(   SWAPRecords(:,3))=   SWAPRecords(:,1);

SWAPIndex=zeros(size(D));

if ~isempty(SWAPRecords)
    % sort SWAP records by tree log L (descending)
    SWAPRecords=sortrows(SWAPRecords,-2);
    % so the first row is the "best" swap, and assign 2 as priority
    SWAPIndex(SWAPRecords(1,3))=2;
    SWAPIndex(SWAPRecords(2:end,3))=1;
end

% build up new CMatrix from SWAPPING records
NewCMatrix=CMatrix;

for i=1:size(SWAPRecords,1)
    
    NodeIndex   =SWAPRecords(i,3);
    ParentIndex =SWAPRecords(i,4);
    ChildIndex  =SWAPRecords(i,5);
    BrotherIndex=SWAPRecords(i,6);

    NewCMatrix(ParentIndex,BrotherIndex)=0;
    NewCMatrix(NodeIndex  ,BrotherIndex)=1;

    NewCMatrix(NodeIndex  ,ChildIndex)=0;
    NewCMatrix(ParentIndex,ChildIndex)=1;
end



function Records=CombineTable(Table_noSWAP,Table_SWAP,TreeLogL)
T=[Table_noSWAP;Table_SWAP];
T(T(:,2)<TreeLogL,:)=[]; % remove nonsense optimization

T=sortrows(T,-2); % debugging use, good for visualize data

% indicator if current node index != parent node index
IsSWAP=T(:,3)~=T(:,4);

Records=cell(size(T,1),1);

for i=1:numel(Records)
    [MaxValue,Index]=max(T(:,2));
    
    % if max value == -inf
    % it means there's no more availiable edge optimization
    % then quick
    if isinf(MaxValue), break; end
    
    % keep records
    Records{i}=T(Index,:);
    
    % this is for SWAP only
    if IsSWAP(Index)
        % remove records with adjacent edge including itself to prevent dependent swapping
        % to guarantee swappings are independent
        Edge=[T(Index,3) T(Index,4)];
        j=IsSWAP&(ismember(T(:,3),Edge)|ismember(T(:,4),Edge));
        T(j,2)=-inf;        
    end
    
    % because node_Index has been picked up, remove all the other potential
    % optimization for node_Index.
    % I use -inf as a trick, since Index is assigned by tree loglikelihood
    j=T(:,3)==T(Index,3);
    T(j,2)=-inf;
    
end

% Records is a condensed table with all possible optimization (swap or non-swap)
Records=vertcat(Records{:});


