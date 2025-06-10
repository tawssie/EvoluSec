function [Dist,LogLikelihood]=AADIST(Model,SeqList)
% This is the exact implementation of PRODIST algorithm
% by Joseph Felsenstein with maximum likelihood approach
% 1. Non-standard amino acid are not considered
% 2. Gap symbol = "-"
% 3. Sequences must be well-aligned, i.e., same length 
%
% % Below are examples from PROML program (cases are different between PRODIST and PROML :~)
% s1='AACGTGGCCAAAT';
% s2='AAGGTCGCCAAAC';
% s3='CATTTCGTCACAA';
% s4='GGTATTTCGGCCT';
% s5='GGGATCTCGGCCC';
% 
% 
% % Below are examples from PRODIST program
% s1='AACGTGGCCACAT';
% s2='AAGGTCGCCACAC';
% s3='CAGTTCGCCACAA';
% s4='GAGATTTCCGCCT';
% s5='GAGATCTCCGCCC';
% 
% SeqList={s1,s2,s3,s4,s5}';
% Model=JTT;
% [Dist,LogLikelihood]=AADIST( Model,{s1,s2,s3,s4,s5});
% DistanceMatrix=squareform(Dist);
%


% "ones" function with classname doen't support old MATLAB version
% [i,j]=ind2sub([numel(SeqList) numel(SeqList)],find(tril(ones(numel(SeqList),'logical'),-1)));
[i,j]=ind2sub([numel(SeqList) numel(SeqList)],find(tril(ones(numel(SeqList))==1,-1)));
Dist         =zeros(numel(i),1);
LogLikelihood=zeros(numel(i),1);
 

D=diag(Model.D);
F=Model.F;
V=Model.V;

for k=1:numel(Dist)
    N=GetObservedChangeMatrix(SeqList{i(k)},SeqList{j(k)});
    [Dist(k), LogLikelihood(k)]=GetEvoDistanceByMLModel(F,V,D,N);
end



function [EvoDist,LogLikelihood]=GetEvoDistanceByMLModel(F,V,D,N)
% minimize "negative" likelihood == maximize likelihood

f=@(t)(-GetLikelihood(F,V,D,N,t));

[d,e]=BrentOptAlgo4Dist(f,0,50);

EvoDist=d;
LogLikelihood=-e;


function N=GetObservedChangeMatrix(s1,s2)
% This function counts AA changes between 2 sequences

isGap=(s1=='-'|s2=='-');
s1(isGap)=[];
s2(isGap)=[];

s1=s1-'A'+1;
s2=s2-'A'+1;

N = histcounts2(s1,s2,1:27,1:27);

N('BJOUXZ'-'A'+1,:)=[];
N(:,'BJOUXZ'-'A'+1)=[];

function L=GetLikelihood(F,V,D,ObservedChangeMatrix,Time)
% P=GetExpMByTime(V,D,time);
% F = normalized frequency
% V = matrix of eigenvectors
% D = eigenvalues
% ObservedChangeMatrix = counts of differences within a pair of sequneces
% Time = evolutionary distance between a pair of sequences
% P = probability matrix calculated by matrix exponential

P=V*diag(exp(D*Time))/V;
L=sum(sum((log(repmat(F,1,20))+log(P)).*ObservedChangeMatrix));

