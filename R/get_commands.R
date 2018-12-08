
##########################################################################
#get parameters from the command line. (e.g. readcount file, bed file ,..)
##########################################################################

args=commandArgs(TRUE)
commands=c()

#store optional parameters in a matrix. if args == 0 do the default settings!
if(length(args)!=0) {
	for(a in 1:length(args)){
		if(length(grep("-(i|o|h)",args[a])) != 0){
			commands=rbind(commands,c(args[a],""))
		}else{
			if(commands[dim(commands)[1],2] == ""){
				commands[dim(commands)[1],2] = args[a]
			}
		}
	}
}else{
	
	#stop if no arguments are used:
	stop(error="\nUsage:\tMDS.R -h -i SimilarityMatrix.txt -o output.pdf \n
		parameters:
		-i\tSimilarity Matrix
		-o\tname of output .pdf file
		
		optional parameters:
		-h\tshows this help message
		\n",call. = FALSE)
}


#stop if the obligatory arguments are not used or the -h option is used:
if(length(which(args == "-i")) == 0 || length(which(args == "-h")) != 0){
	
	#stop if the obligatory arguments are not used:
	stop(error="\nUsage:\tMDS.R -h -i SimilarityMatrix.txt -o output.pdf \n	
		parameters:
		-i\tSimilarity Matrix
		-o\tname of output .pdf file
		
		optional parameters:
		-h\tshows this help message
		\n",call. = FALSE)
}

#get similarity matrix:
test_sample=commands[which(commands[,1] == "-i"),2]

#get name for output pdf file:
output=commands[which(commands[,1] == "-o"),2]
