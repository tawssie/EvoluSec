function [Dist,LogLikelihood]=SSDIST(Model,SeqList)
% Pairwise evolutionary distance estimation by maximum likelihood approach
% with secondary structure model.
% 1. Only accepts "BEGHIST_-" labels
% 2. "BEGHIST" are defined by DSSP program.
% 3. The "_" label represents "unknown" secondary structure, it is NOT gap
% 3. Gap symbol = "-"
% 4. Secondary structure must be well-aligned, i.e., same length
% "UNKNOWN secondarey structure status" are treated in a special way
%     because "UNKNOWN status" is NOT "Gap"
%
% eg., example 1 with few unknown status
% s1='GHBETGIBSHBHT__III__';
% s2='H_GGTBGBBHBHBEE_____';
% s3='BHGTTBGBBHBHHBGBBGBI';
% s4='GHGHTTTBBGBBTBGBBGBI';
% s5='GHGHTBTBBGBBBBGBBGBI';

% eg., example 2 secondary structure are all known
% s1='GBETGIBSHBHT';
% s2='HGGTBGBBHBHB';
% s3='BGTTBGBBHBHH';
% s4='GGHTTTBBGBBT';
% s5='GGHTBTBBGBBB';

% SeqList={s1,s2,s3,s4,s5}';
% Model=EvoSS('DSO-SS');
%
% [Dist,LogLikelihood]=SSDIST(Model,{s1,s2,s3,s4,s5});
%
% DistanceMatrix=squareform(Dist);
%

% "ones" function with classname doen't support previous MATLAB version
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
 
[d,e]=BrentOptAlgo4Dist(f,0,50); % customised optimisation

EvoDist=d;
LogLikelihood=-e;


function N=GetObservedChangeMatrix(s1,s2)
isGap=(s1=='-'|s2=='-');
s1(isGap)=[];
s2(isGap)=[];

s1=s1-'A'+1;
s2=s2-'A'+1;
                           
SSOrder='BEGHIST_'-'A'+1;
N_BEGHIST_=histcounts2(double(s1),double(s2),[double(SSOrder) inf],[double(SSOrder) inf]);

N=N_BEGHIST_(1:7,1:7); % N records changes within BEGHIST 7 states


N_Unknown1=repmat(N_BEGHIST_(8,1:7),[7,1])/7;
N_Unknown2=repmat(N_BEGHIST_(1:7,8),[1,7])/7;
N_Unknown3=repmat(N_BEGHIST_(8,8),[7,7])/(7*7);

N=N+N_Unknown1+N_Unknown2+N_Unknown3;


function L=GetLikelihood(F,V,D,ObservedChangeMatrix,Time)
% F = normalized frequency
% V = matrix of eigenvectors
% D = eigenvalues
% ObservedChangeMatrix = counts of differences within a pair of sequneces
% Time = evolutionary distance between a pair of sequences
% P = probability matrix calculated by matrix exponential

P=V*diag(exp(D*Time))/V;
L=sum(sum((log(repmat(F,1,7))+log(P)).*ObservedChangeMatrix));

