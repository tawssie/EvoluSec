function EstimateSSDistance(InputFileName,OutputFileName)
% USAGE:
% This function estimates evolutionary distance based on
% secondary-structure evolutionary model (DSO-SS).
% Input file is multiple secondary structure alignment in fasta format.
%
% E.g., InputFileName ='example_PF00017_SS.fasta';
%       OutputFileName='example_PF00017_SS_result.txt';
%
% EstimateSSDistance('example_PF00017_SS.fasta','example_PF00017_SS_result.txt');


Main;

Data=fastaread(InputFileName);

Model=EvoSS('DSO-SS');
D=SSDIST(Model,{Data.Sequence});

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


