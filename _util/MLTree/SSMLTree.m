function [D,TreeLogL,CMatrix,ID,ExtraResult]=SSMLTree(MSSA,InitTreeAlgo,Model,TreeLogLStableThreshold,DoParfor)
% function [D,TreeLogL,CMatrix,ID,ExtraResult]=SSMLTree(MSSA,InitTreeAlgo,DoParfor)

% This is the PHYML version to optimize Tree likelihood
% by improving all edges with lambda 0.75 iteratively
% and subtree SWAPPING tree procedure.
% This function uses customized fminbnd to solve extreme values.
%
% Secondary structure version
%
% Use Adachis' least square calculation to assign initial edge (branch) lengths of topology
% Use JF's conditional likelihood calculcation (prunning)
% Iteratively optimize all edges until tree likelihood stable with
% threshold


if ~(nargin==4||nargin==5)
    error('not enough arguments, see the usage');
end


OutputFile=fopen('RunningLog.txt','w');


ExtraResult=cell(400,1); % Used to store tree information when swap happens
EResult.CMatrix=[];      % Topology
EResult.D=[];            % Edge length attached to each node
EResult.TreeLogL=[];     % Log likelihood for whole tree
EResult.SWAPIndex=[];    % Swaps happening at the node (index)


% 1. Choose evolutionary model
% Model=EvoSS('DSO-SS'); 

% 2. Build up tree by NJ or UPGMA with estimated (ML) pairwise distance
Tree=BuildInitTree(Model,MSSA,InitTreeAlgo);

% 3. Retrieve tree information:
%    connected matrix, root index, traverse order (reverse of DFS)
[CMatrix,ID]=getmatrix(Tree);
CMatrix=full(CMatrix);
RootIndex=find(~any(CMatrix,1));% An easy way is to find a column with no inticator

TraverseOrder=graphtraverse(sparse(CMatrix),RootIndex,'method','DFS');
TraverseOrder=flipud(TraverseOrder(:));



% 4. Make sure, MSA header ordering == tips' name ordering
%     ***IMPORTANT IMPORTANT IMPORTANT***
%     MSA order must be consistent with tree tips' name orderering
%     i.e., MSA order must follow the first N tips ID
%     because observed conditional vector for tips 
%     is assigned according to MSA's order
[~,i]=ismember(ID(1:numel(MSSA)),{MSSA.Header});
MSSA=MSSA(i);clear i;



% 5. Use constraint least square method to intialize edge length.
%     Since order of MSA has been changed according to tips' node name
%     , there's a need to rearrange pairwise distance.
PairwiseDist=SSDIST(Model,{MSSA.Sequence});% by ML estimation

D_initial=GetInitialDistanceByLeastSquare(CMatrix,PairwiseDist);

D_initial(D_initial<1e-6)=1e-6;
D_initial(end)=[]; % Since the last element is
                   % "pseudo" edge length of root to its parent
                   % but root does't have its parent, so remove it
                   % ps. D_initial(5) == edge length connects to node 5
D=D_initial;


% 6. Transform sequence alignment (ARNK..) to MSALabel (matrix of number)
%    Define TipLikelihood by Model
[TipLikelihood,MSSALabel]=GetTipsLikelihood4SS(vertcat(MSSA.Sequence),Model);

% 7. Calculate initial tree log likelihood
%    Preserve all conditional likelihood for all nodes, all children
%    Preserve all transition probability matrix
[TreeLogL,CondL,P]=GetTreeLikelihood03(TipLikelihood,MSSALabel,Model,CMatrix,TraverseOrder,D);


% Pre calculate lambda parameters, starting by 0.75 and divided by 2, 0.75*2^-10 at most
LambdaList=[0.75*power(2,0:-1:-10) 0]; % miminum NON-ZERO lambda=0.000732421875

IterationTree=0;           % # of iter to update a tree (remember topology and distances are done altogether)
IterationLambda=0;         % # of iter to use smaller lambda, i.e., lambda != 0.75
MaxIterationTree=200;      % To control maximum iteration
% TreeLogLStableThreshold=1e-8;


while true

    % Display current tree information
    fprintf('Tree Log L =%.5f, Tree iteration=%05d, Lambda iteration=%05d\n',TreeLogL,IterationTree,IterationLambda);
    fprintf(OutputFile,'Tree Log L =%.5f, Tree iteration=%05d, Lambda iteration=%05d\n',TreeLogL,IterationTree,IterationLambda);

    % Retrieve new
    % 1) Topology 2) Potential better edge lengths 3) Traverse order
    [CMatrix,NewD,SWAPIndex]=ImproveTree04(TreeLogL,CondL,P,Model,CMatrix,D,DoParfor);
    
    % ImproveTree will provide new topology, so TraverseOrder must be updated
    TraverseOrder=graphtraverse(sparse(CMatrix),RootIndex,'method','DFS');
    TraverseOrder=flipud(TraverseOrder(:));

    % Go through all lambda untill a better NewTreeL or reach the end of LambdaList
    % (Lambda=Zero) means D won't get any change from NewD
    for Lambda=LambdaList
        
        % Update distance by Lambda
        LambdaNewD=LambdaModification(Lambda,D,NewD,SWAPIndex);
        
        if any(LambdaNewD<1e-6)
            LambdaNewD(LambdaNewD<1e-6)=1e-6;
        end
        
        % Estimate New Tree Log L from new distance
        [NewTreeLogL,CondL,P]=GetTreeLikelihood03(TipLikelihood,MSSALabel,Model,CMatrix,TraverseOrder,LambdaNewD);
        
        % If New Tree Log L is better than previous Tree Log L than leave
        % Otherwise, try different Lambda
        if (NewTreeLogL>=TreeLogL)
            break;
        end
        
        IterationLambda=IterationLambda+1;
    end
    
    % It means Lambda hits zero, may not find a better lambda
    if Lambda==0
        IsEndless=true;
        warning('Lambda reaches pre-defined smallest value as %.14f',LambdaList(end-1));
    else
        IsEndless=false;
    end
    
    % Test whether New Tree L is stable by threshold
    IsStable=abs(NewTreeLogL-TreeLogL)<TreeLogLStableThreshold;
    
    % Test whether New Tree Log L is better previous Tree Log L
    IsBetterLogL=NewTreeLogL>=TreeLogL;
    
    if IsBetterLogL
        D=LambdaNewD;
        TreeLogL=NewTreeLogL;
    end
       
    IterationTree=IterationTree+1;
    
    % If swap happens, store it to ExtraResult
    if ~all(SWAPIndex==0)
        EResult.CMatrix=CMatrix;
        EResult.D=D;
        EResult.TreeLogL=TreeLogL;
        EResult.SWAPIndex=SWAPIndex;
        ExtraResult{IterationTree}=EResult;
    end
    
    % If Tree Log L is stable, then terminate procedure
    if IsStable
        break;
    end
    
    if IterationTree==MaxIterationTree
        warning('Terminated, reach MaxIterationTree=%03d',IterationTree);
        break;
    end
    
    if ~IsStable&&IsEndless&&~IsBetterLogL
        warning('***Caution*** cannot find a better lamba to improve tree edges at iteration %03d'...
            ,IterationTree);
        break;
    end
end

ExtraResult(IterationTree+1:end)=[];

D=[D;0];

fprintf('Total tree iteration=%05d\nTotal Lambda iteration=%05d\n',IterationTree,IterationLambda);
fprintf('Current Output Tree Log Likelihood=%.8f\n\n',NewTreeLogL);

fprintf(OutputFile,'Total tree iteration=%05d\nTotal Lambda iteration=%05d\n',IterationTree,IterationLambda);
fprintf(OutputFile,'Current Output Tree Log Likelihood=%.8f\n',NewTreeLogL);
fprintf(OutputFile,'Initial Tree is built by %s\n',InitTreeAlgo);
fprintf(OutputFile,'Tree Log Likelihood Threshold =%.8f\n',TreeLogLStableThreshold);

fclose(OutputFile);

 

function LambdaNewD=LambdaModification(Lambda,D,NewD,SWAPIndex)
% SWAPIndex == 2 best SWAP, no    lambda applied
% SWAPIndex == 1 good SWAP, apply lambda proportion
% SWAPIndex == 0 no   SWAP, apply l+lambda*(Newl-l)
% LambdaNewD=NewD;
LambdaNewD=zeros(size(NewD));

i=SWAPIndex == 2; %%% no changes from NewD
LambdaNewD(i)=NewD(i);

i=SWAPIndex == 1;

LambdaNewD(i)= Lambda*(NewD(i)-D(i))+D(i);

i=SWAPIndex == 0;
LambdaNewD(i)= Lambda*(NewD(i)-D(i))+D(i);


function [Tree,PairwiseDist]=BuildInitTree(Model,MSSA,InitTreeAlgo)
% 1. Calculate pairwise distance for protein sequences
PairwiseDist=SSDIST(Model,{MSSA.Sequence});% by ML estimation

if strcmpi(InitTreeAlgo,'NJ')
    % 2.1 Build up initial tree by NJ
    Tree=seqneighjoin(PairwiseDist,'equivar',{MSSA.Header});
elseif strcmpi(InitTreeAlgo,'UPGMA')
    % 2.2 Or build up initial tree by UPGMA
    Tree = seqlinkage(PairwiseDist,'UPGMA',{MSSA.Header});
end

