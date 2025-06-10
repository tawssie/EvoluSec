function EstimateAADistance(InputFileName,OutputFileName)
% USAGE:
% This function estimates evolutionary distance based on
% amino-acid evolutionary model (JTT).
% Input file is multiple sequence alignment in fasta format.
%
% E.g., InputFileName ='example_PF00017_AA.fasta';
%       OutputFileName='example_PF00017_AA_result.txt';
%
% EstimateAADistance('example_PF00017_AA.fasta','example_PF00017_AA_result.txt');


Main;

Data=fastaread(InputFileName);

Model=JTT;

D=AADIST(Model,{Data.Sequence});

D=squareform(D);


File=fopen(OutputFileName,'w');

for i=1:numel(Data)
    fprintf(File,'\t%s',Data(i).Header);
end
fprintf(File,'\n');

NCol=size(D,2);

for i=1:numel(Data)
    fprintf(File,'%s',Data(i).Header);
    
    for j=1:NCol
        fprintf(File,'\t%.8f',D(i,j));
    end
    fprintf(File,'\n');
end
fclose(File);


