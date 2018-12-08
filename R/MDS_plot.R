
##########################################################################################################################################################################
#Create a full Similarity Matrix:
 	Matrix[,(dim(Genomes)[1]+1):(dim(Genomes)[1]+number_samples)]<-t(Matrix[(dim(Genomes)[1]+1):(dim(Genomes)[1]+ number_samples),])
 	Matrix[1:dim(Genomes)[1],1:dim(Genomes)[1]]<-Genomes       	
 	Matrix<-as.matrix(Matrix)
 	
 	INFO_tmp = scan(	test_sample, 
					nlines=1,
					what=character(),
					quiet=T,
					sep="\n",
					skip=1
				)
				
  	write.table(INFO_tmp,file=sub(".pdf","_full_matrix.txt",output),quote=F,sep="\t",col.names=F,row.names=F)
 	write.table(Matrix,file=sub(".pdf","_full_matrix.txt",output),quote=F,sep="\t",col.names=F,row.names=F,append=T)

##########################################################################################################################################################################
#Non-Metric Multidimensional Scaling with isoMDS{MASS}
 	coord<-isoMDS(Matrix,y=cmdscale(Matrix),trace=F)  
  
	##########################################################################################################################################################################
	#Determine zoom factor for plot::
	zoom=0.3
	dev_x=(max(coord$points[,1])-min(coord$points[,1]))*zoom		
	dev_y=(max(coord$points[,2])-min(coord$points[,2]))*zoom

##########################################################################################################################################################################
#Create an empty plot:
	emptyplot(	xlim = c(min(coord$points[,1])-dev_x, max(coord$points[,1])+dev_x), 
  				ylim = c(min(coord$points[,2])-dev_y, max(coord$points[,2])+dev_y),
  				col=gray(0.99),
  				asp=(max(coord$points[,1])-min(coord$points[,1]))/(max(coord$points[,2])-min(coord$points[,2]))
  			)	
	
	#title and axes:
	title(main="Distance to high quality reference set",cex.main=2.5,col.main="black",line=1)
	mtext(expression(phi[1]),side=1,cex=3.0,col="black",line=1.5)
  	mtext(expression(phi[2]),side=2,cex=3.0,col="black",line=0)
  	
  	#add points:
		#1000 Genomes Individuals:
		points(	coord$points[1:dim(Genomes)[1],1],
				coord$points[1:dim(Genomes)[1],2],
				col="black",
				pch=1,
				cex=4
			)

		#Test Samples:
		for(i in 1:number_samples){
  			points(	coord$points[(dim(Genomes)[1]+i),1],
  					coord$points[(dim(Genomes)[1]+i),2],
  					col= "firebrick",
  					pch=2,
  					cex=5
  				)
  			text(	coord$points[(dim(Genomes)[1]+i),1],
  					coord$points[(dim(Genomes)[1]+i),2],
  					label=(i),
  					col="firebrick",
  					pch=2,
  					cex=2
  				)		
  		}
  	
  	#add inner legend:
  	legend(	x=min(coord$points[,1])-dev_x,
  			y=max(coord$points[,2])+dev_y,
  			legend=c(paste("1000Genomes ",paste(population,"*",sep="")," Individuals",sep=""),
  			paste("Test Samples",sep="")),
  			col=c("black","firebrick"),
  			cex=2,
  			bty="n",
  			text.col=c("black","firebrick"),
  			pch=c(1,rep(2, number_samples))
  		)
	
	#population description:
	pop_descr=Ref[index+10]
	legend(	"bottomleft",
			legend=paste("*",pop_descr,sep=""),
  			cex=1.3,
  			bty="n",
  		)
