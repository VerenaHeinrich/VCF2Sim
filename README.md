
# *********************  VCF2Sim  ********************** #


##########################################################
# get the similarity values for all samples in a VCF file:
##########################################################

Run 'java -jar VCF2Sim.jar' to the all parameters:

usage: VCF2Sim.jar
 -d,--dir <arg>      Directory to vcf files
 -g,--gfv <arg>      Name of optional GFV file. default build on 1KGP
                     ('/VCF2Sim/data/GFV_consensus.txt')
 -m,--metric <arg>   Metric (Sensitive for 'rare' or 'common' variants.
                     You can also use a normal 'hamming' distance. default
                     = 'rare')
 -o,--out <arg>      Name of output file
 -p,--pop <arg>      1000Genomes background population (if not stated, the
                     most likely will be selected.)
                     

Run 'java -jar VCF2Sim.jar  -d <folder to .vcf files> -o <output file>'
to get similarities with weight on rare variants. 

Output files are in the input directory.


##########################################################
# visualize the results and get the quality report:
##########################################################

Run 'Rscript MDS.R' to the all parameters:

Usage:	MDS.R -h -i SimilarityMatrix.txt -o output.pdf 
	
		parameters:
		-i	Similarity Matrix
		-o	name of output .pdf file
		
		optional parameters:
		-h	shows this help message



