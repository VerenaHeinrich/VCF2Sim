############################
#Similarity - Visualization#
############################

####################################################################################
#packages:
	library(MASS)   #for MDS
	library(shape)  #for creating an emptyplot
	
#################################################################################### 
#get commands:
	source("R/get_commands.R")
	
	#background population (e.g.: CEU)
	INFO = scan(	test_sample, 
					nlines=1,
					what=character(),
					quiet=T,
					sep="\t",
					skip = 1
				)	
	population  =""

####################################################################################
#Input Data:
 	Matrix<-read.table(test_sample,header =FALSE)

  	if(grep("^#",INFO[1])==1){
    	population = sub("#","",INFO[1])
 	}else{
	  	if(length(args == 3)){population=args[3]}
	 	 else(print("Header is missing! Please specify a background population!"))
  			stop("\nUsage:\tMDS.R dir/test_sample.txt dir/output\n\n")
	}
  	
  	INFO=as.data.frame(strsplit(INFO[2:length(INFO)],":"),stringsAsFactors=F)
  	for(i in 1:length(INFO)){
  		INFO[[i]][1] =(gsub("_"," ",INFO[[i]][1]))

  		if(nchar(as.character(INFO[[i]][1]))>10){
  			INFO[[i]][1]=paste(substr(INFO[[i]][1],1,10),"[..]",sep="")
  		}
  		for(j in 2:6){
  			INFO[[i]][j] = round(as.numeric(as.character(INFO[[i]][j])),digits=3)
  			if(as.character(INFO[[i]][j]) == "NaN"){
  				INFO[[i]][j] = "-"
  			}
  			else if(nchar(as.character(INFO[[i]][j]))<5){
  				INFO[[i]][j]=paste(INFO[[i]][j],"0",sep="")
  			}
  		}
  	}
  	
  	Genomes<-read.table(	paste("Matrices/D_G3_1000genomes_consensus_", 
  							population,".txt",sep=""),
  							header=FALSE
  						)
 	number_samples<-dim(Matrix)[1]-dim(Genomes)[1]
  
####################################################################################
#Read Reference file:
Ref<-scan("Reference_Curve.txt",what=character(),sep="\n",quiet=T)
index<-which(Ref==paste("#",population,sep=""))

####################################################################################
#MDS plot:
	pdf(sub(".pdf","_MDS.pdf",output), height =10,width=10)
			par(oma=c(0,0,0,0),mar=c(5,5,5,5))
		  	source("R/MDS_plot.R")
	dev.off()
	
####################################################################################
#SDS plot:
	pdf(sub(".pdf","_SDS.pdf",output), height =10,width=10)
			par(oma=c(0,0,0,0),mar=c(5,5,5,5))
		  	source("R/SDS_plot.R")
	dev.off()
	
####################################################################################
#Read Latex template:
	tex=scan(file="LATEX/report.tex",sep="\n",what="character")

####################################################################################
#change template file:
	tex=sub("ValueForDirectory",sub("/.*$","",output),tex)
	tex=sub("ValueForNumberOfSamples",number_samples,tex)
	
	#value MDS plot:
	tex=sub("ValueForMDSPlot",sub(".pdf","_MDS.pdf",output),tex)

	#value for SDS plot:
	tex=sub("ValueForSDSPlot",sub(".pdf","_SDS.pdf",output),tex)

	#sample names:
	ValueForSamples=""
	for(i in 1:length(unlist(INFO[1,]))){
		ValueForSamples=paste(ValueForSamples,"(",i,") ",as.character(unlist(INFO[1,i])),sep="")
		if(i<length(unlist(INFO[1,]))){
			ValueForSamples=paste(ValueForSamples,"&",sep="")
		}
	}
	tex=sub("ValueForSamples",as.character(ValueForSamples),tex)
	
	#values for ti/tv-ratio
	ValueForTrans=""
	for(i in 1:length(unlist(INFO[1,]))){
		ValueForTrans=paste(ValueForTrans,INFO[1:dim(INFO)[2]][2,i], sep="")
		if(i<length(unlist(INFO[1,]))){
			ValueForTrans=paste(ValueForTrans,"&",sep="")
		}
	}
	tex=sub("ValueForTrans",ValueForTrans,tex)

	#values for dbsnp
	ValueForDBSNP=""
	for(i in 1:length(unlist(INFO[1,]))){
		ValueForDBSNP=paste(ValueForDBSNP,INFO[1:dim(INFO)[2]][3,i], sep="")
		if(i<length(unlist(INFO[1,]))){
			ValueForDBSNP=paste(ValueForDBSNP,"&",sep="")
		}
	}
	tex=sub("ValueForDBSNP",ValueForDBSNP,tex)

	#values for var(het)
	ValueForVariance=""
	for(i in 1:length(unlist(INFO[1,]))){
		ValueForVariance=paste(ValueForVariance,INFO[1:dim(INFO)[2]][4,i], sep="")
		if(i<length(unlist(INFO[1,]))){
			ValueForVariance=paste(ValueForVariance,"&",sep="")
		}
	}
	tex=sub("ValueForVariance",ValueForVariance,tex)
	
	#values for SDS
	ValueForSDS=""
	for(i in 1:length(unlist(INFO[1,]))){
		ValueForSDS=paste(ValueForSDS,round(median[i],digits=3), sep="")
		if(i<length(unlist(INFO[1,]))){
			ValueForSDS=paste(ValueForSDS,"&",sep="")
		}
	}
	tex=sub("ValueForSDS",ValueForSDS,tex)

	#values for Nonsynonymous/Synonymus Ratio:
	ValueForSynTrans=""
	for(i in 1:length(unlist(INFO[1,]))){
		ValueForSynTrans=paste(ValueForSynTrans,INFO[1:dim(INFO)[2]][5,i], sep="")
		if(i<length(unlist(INFO[1,]))){
			ValueForSynTrans=paste(ValueForSynTrans,"&",sep="")
		}
	}
	tex=sub("ValueForSynTrans",ValueForSynTrans,tex)

	#values for het/hom ratio:
	ValueForHetHomRatio=""
	for(i in 1:length(unlist(INFO[1,]))){
		ValueForHetHomRatio=paste(ValueForHetHomRatio,INFO[1:dim(INFO)[2]][6,i], sep="")
		if(i<length(unlist(INFO[1,]))){
			ValueForHetHomRatio=paste(ValueForHetHomRatio,"&",sep="")
		}
	}
	tex=sub("ValueForHetHomRatio",ValueForHetHomRatio,tex)

####################################################################################
#Replace template with current parameters and write out:
	write(tex, file=sub(".pdf",".tex",output))
	
####################################################################################
#Change directory and execute Latex .tex file:
	
	tmp=strsplit(output,"/")	
	out_dir=""
	for(i in 1:(length(unlist(tmp))-1)){
		out_dir=paste(out_dir,unlist(tmp)[i],"/",sep="")
	}
	
		
	table_out=c()
	for(t in 1:length(unlist(INFO[1,]))){	
		table_out=cbind(table_out,c(	INFO[1:dim(INFO)[2]][2,t],
										INFO[1:dim(INFO)[2]][3,t],
										INFO[1:dim(INFO)[2]][4,t],
										INFO[1:dim(INFO)[2]][5,t],
										INFO[1:dim(INFO)[2]][6,t],
										round(median[t],digits=3)
									)
						)
		
	}
	rownames(table_out) =c(	"ti/tv-ratio",
							"Percentage in dbSNP138",
							"var(het)",
							"Nonsynonymous/Synonymus Ratio",
							"het/hom ratio",
							"SDS")
	colnames(table_out) = unlist(INFO[1:dim(INFO)[2]][1,])
 	write.table(	(table_out),
					file=sub(".pdf","_table.txt",output),
					sep="\t",
					row.names=T,
					col.names=T,
					quote=F
				)
	
	system(paste("pdflatex --shell-escape --interaction=nonstopmode --interaction=batchmode -aux-directory=\"",out_dir,"\" -output-directory=\"",out_dir,"\" \"",sub(".pdf",".tex",output),"\"",sep=""))
	