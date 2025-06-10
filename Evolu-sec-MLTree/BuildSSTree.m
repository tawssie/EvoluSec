function BuildSSTree(InputMSSAFileName,DoParfor)
% This function provides simple procedure to construct a phylogenetic tree 
%     based on Evolu-sec and multiple secondary structure alignment.
% USAGE:
%     BuildSSTree(InputMSSAFileName, DoParfor);
% Input:
%     InputMSSAFileName 
%         is a fasta file that contains multiple secondary structure alignment.
%     DoParfor
%         is a parameter to control parallel computing
%         (needs Parallel Computing Toolbox)
% Output:
%     Two files are generated for the secondary structure phylogenetic tree.
%     They are InputMSSAFileName".tree" and InputMSSAFileName".log".
%     eg. InputMSSAFileName="example.fasta"
%     "example.fasta.tree" records the inferred tree as Newick string.
%     "example.fasta.log"  records the log likelihood for each iteration.
% 
% Example:
%     BuildSSTree('example.fasta');       % default is without parallel computing
%     BuildSSTree('example.fasta',false); % without parallel computing
%     BuildSSTree('example.fasta',true);  % use parallel computing



Main;
if ~exist('DoParfor','var')
    DoParfor=false;
end

MSSA=fastaread(InputMSSAFileName);

[D,TreeLogL,CMatrix,ID]=SSMLTree(MSSA,'UPGMA',EvoSS('DSO-SS'),1e-8,DoParfor);

PhyTreeObj=CreatePhytreeFromBGraph(CMatrix,ID,D);
NWString=getnewickstr(PhyTreeObj);

File=fopen([InputMSSAFileName '.tree'],'w');
fprintf(File,'%s\n',NWString);
fclose(File);

movefile('RunningLog.txt',[InputMSSAFileName '.log']);



