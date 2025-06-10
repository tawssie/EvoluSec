USAGE:
01. Change MATLAB working directory to "Evolu-sec-Package/Evolu-sec-MLTree"
02. Use "Main" in Matlab terminal to add paths for supporting scripts/functions
03. Type "help BuildSSTree" in Matlab terminal to show additional information
04. Type following 1 command in MATLAB to run an example for secondary structure phylogenetic tree inference
    "BuildSSTree('example.fasta');" % note ***
    
05. Results are calculated and exported in
    "example.fasta.tree"
    and
    "example.fasta.log"

% note ***
	Default is without parallel computing.
	It may need "MATLAB Parallel computing toolbox" to speed up the inference.
	Old MATLAB version for parallel computing may need additional requirement to activate parallel computation.	
 