function DoASSR(InputMSSAFileName, InputFilePhyTree)
% This function provides simple interface to perform 
%     ancestral secondary structure reconstruction.
% USAGE:
%     DoASSR(InputMSSAFileName, InputPhyTreeFileName);
% Input:
%     InputMSSAFileName 
%         is a fasta file that contains multiple secondary structure alignment.
%     InputPhyTreeFileName
%         is a text file that describe phylogenetic tree.
%         This tree could be NEXUS or NEWICK format.
%         Leaves nodes must contain names that identitcal to header of MSSA.
%         Branch(internal) nodes must provide names to specificy ancestors.
% Output:
%     Two output files are generated for ASSR predictions.
%     They are InputMSSAFileName".ancSS.txt" and InputMSSAFileName".log"
%
%     eg.
%     InputMSSAFileName='example_Wilson_MSSA.fasta';
%     InputFilePhyTree='example_Wilson_PrunedTree.txt';
%     DoASSR(InputMSSAFileName, InputFilePhyTree);
%
%     "example_Wilson_MSSA.fasta.ancSS.txt" contains inferred secondary structure.
%     "example_Wilson_MSSA.fasta.prob.txt" contains joint probabilities.


Main;


OutputSSFileName=[InputMSSAFileName '.ancSS.txt'];
OutputProbFileName=[InputMSSAFileName '.prob.txt'];

MSSA=fastaread(InputMSSAFileName);
T=phytreeread(InputFilePhyTree);

[MSSA.('SS')] = MSSA.('Sequence');
MSSA = rmfield(MSSA,'Sequence');

[AllCGTableJoint,AllNodeResult,ID]=ASSR(MSSA,T);

SS=cellstr(AllNodeResult);
File=fopen(OutputSSFileName,'w');

for i=1:numel(ID)
    fprintf(File,'>%s\n%s\n',ID{i},SS{i});
end

fclose(File);


SS=cellstr(AllNodeResult');
File=fopen(OutputProbFileName,'w');

for i=1:numel(SS)
    fprintf(File,'%s\t%.8f\n',SS{i},AllCGTableJoint(i));
end

fclose(File);
