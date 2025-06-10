function [AllCGTableJoint,AllNodeResult,ID]=ASSR(MSSA,PhyloTree)
% The function performs ancestral "secondary structure" reconstruction.
% It takes a phylogenetic tree object 
%    and a multiple secondary structure alignment.
% USAGE:
% [AllCGTableJoint,AllNodeResult,ID]=ASR_SS03(MSSA,PhytreeObj);
%
% Fieldnames in MSSA
%     'Header'    = Provided by authors
%     'SS'        = secondary structure labels of B,E,G,H,I,S,T,_ or - (the gap)

[ConnectedMatrix, ID, Dist] = getmatrix(PhyloTree);

% NodeList is created by order of "getmatrix",
% i.e., 1:numel(MSSA) are leaves

% Below is to check consistency between phytree ID and MSSA header
[~,i]=ismember(ID(1:numel(MSSA)),{MSSA.Header});

MSSA=MSSA(i);

if ~isequal({MSSA.Header}',ID(1:numel(MSSA)))
    warn('not equal order');
end

ConnectedMatrix=full(ConnectedMatrix); % use full version first

TimeMatrix=zeros(size(ConnectedMatrix));
TimeMatrix(ConnectedMatrix)=Dist(1:end-1);

LeafInstance=vertcat(MSSA.SS);


NodeList=CreateInferTree(ConnectedMatrix,TimeMatrix,ID);

AllCGTableJoint=cell(size(LeafInstance,2),1);
AllNodeResult=cell(size(LeafInstance,2),1);

for i=1:size(LeafInstance,2)

    if ~any(LeafInstance(:,i)=='-')
        [AllCGTableJoint{i},AllNodeResult{i}]=InferAllCGTable(NodeList,LeafInstance(:,i));
    else
        AllCGTableJoint{i}=[];
        AllNodeResult{i}=repmat('-',size(NodeList));
    end
    
end

AllNodeResult=horzcat(AllNodeResult{:});
AllCGTableJoint=AllCGTableJoint2Double(AllCGTableJoint);


function JointProb=AllCGTableJoint2Double(AllCGTableJoint)

StringList=cellfun(@char,AllCGTableJoint,'uni',0);

[~,StringList]=strtok(StringList,'=');
JointProb=str2double(strtok(StringList,'='));


function [NodeList,VariableList]=CreateInferTree(ConnectedMatrix,TimeMatrix,ID)
% Import all package from bnkit
import bn.node.*;
import bn.*;
import dat.*;
import bn.prob.*;
import bn.alg.*;
import bn.ctmc.*;


% It needs CharVector2JavaArray function to use java Char array
Label='BEGHIST'; % secondary structure label
Label=CharVector2JavaArray(Label);

Model=EvoSS('DSO-SS');


VariableList=cell(size(ID));
for i=1:numel(VariableList)
    VariableList{i}=EnumVariable(Enumerable(Label),ID{i});        
end
VariableList=vertcat(VariableList{:});


NodeList=cell(size(ID));

% Define all probabilities and dependencies
for i=1:size(ConnectedMatrix,2)
    
    ParentList=ConnectedMatrix(:,i);    
    if any(ParentList)

        % Create nodes' dependency by CPT
        NodeList{i}=CPT(VariableList(i),[VariableList(ParentList(:)') '']);
        
        Time=TimeMatrix(ParentList,i);
        
        % GetExpMByTime to get probability by time
        Prob=GetExpMByTime(Model.V,Model.D,Time);
        
        % SetTransitionProb
        NodeList{i}=SetTransitionProb(false,NodeList{i},Enumerable(Label),Prob);     
    else
        % no childlist == prior        
        NodeList{i}=CPT(VariableList(i));
        NodeList{i}=SetTransitionProb(true ,NodeList{i},Enumerable(Label),Model.F);
    end
end

NodeList=vertcat(NodeList{:});


function [CGTableJoint,NodeResult]=InferAllCGTable(NodeList,LeafInstance)
% Import all package from bnkit
import bn.node.*;
import bn.*;
import dat.*;
import bn.prob.*;
import bn.alg.*;
import bn.ctmc.*;

% CAUTION: All nodes MUST clear instance before "inference and reuse"
for i=1:numel(NodeList)
    NodeList(i).resetInstance();
end


% NodeList is created by order of "getmatrix",
% i.e., 1:numel(MSSA) are leaves
for i=1:numel(LeafInstance)
    if LeafInstance(i)=='_'
        continue;
    end
    NodeList(i).setInstance(LeafInstance(i));
end

BN=BNet();
BN.add(NodeList);

VE = VarElim();
VE.instantiate(BN);

QueryJoint = VE.makeMPEforAll();

CGTableJoint = VE.infer(QueryJoint);
ResultAssignmentArray = CGTableJoint.getMPE();

for i=1:numel(ResultAssignmentArray)
    ResultVariable = ResultAssignmentArray(i).var;
    BN.getNode(ResultVariable).setInstance(ResultAssignmentArray(i).val);            
end

NodeResult=arrayfun(@getInstance,NodeList);






