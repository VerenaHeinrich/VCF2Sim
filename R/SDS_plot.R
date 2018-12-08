##########################################################################################################################################################################
#Standardized Disimilarity Score(SDS): Comparison with Reference Curve:
median<-c()
  	for(l in (dim(Genomes)[1]+1):(dim(Matrix)[1])){    
  	 	median<-c(median, median(Matrix[c(1:dim(Genomes)[1],l),l],na.rm=T))
  	}
  	
  	Error_rate<-c("0","0.00001","0.0001","0.001","0.01")
  	Ref<-scan("Reference_Curve.txt",what=character(),sep="\n",quiet=T)
  	index<-which(Ref==paste("#",population,sep=""))
  	ref_median=as.numeric(unlist(Ref[index+1]))
  	ref_IQR=as.numeric(unlist(Ref[index+2]))
  
 	median<-(median-ref_median)/sqrt(ref_IQR)		#(median-reference_median)/sqrt(reference_IQR)
  	Zscore<-as.numeric(unlist(strsplit(Ref[index+3]," ")))
  	konf<-as.numeric(unlist(strsplit(Ref[index+4]," ")))
  	colors=rainbow(start=0.1,end=1,length(median))


##########################################################################################################################################################################
#Create an empty plot:
	
	emptyplot(	xlim = c(1-0.1,length(Error_rate)+0.1),
    	        ylim = c(-0.1, max(Zscore)+0.1),
           		col=gray(0.99),
           		asp=(length(Error_rate)+0.1-1+0.1)/((max(Zscore))+0.1+0.1)
            )
            
   	#title and axes:
  	title(	main="Genotyping Accuracy",
  			cex.main=2.5,
  			col.main="black",
  			line=1
  		)
	axis(side=1,at=seq(1:length(Error_rate)),labels=(as.character(Error_rate)),cex.axis=2,line=NA,col="dimgray",col.ticks="dimgray",lwd=0,lwd.ticks=1)
  	axis(side=2,at=seq(0,(max(Zscore)+0.1),by=0.5),labels=seq(0,(max(Zscore)+0.1),by=0.5),cex.axis=2,line=NA,col="dimgray",col.ticks="dimgray",lwd=0,lwd.ticks=1)

  	mtext("Estimated Error Rates",side=1,cex=2,col="black",line=3.5)
  	mtext("Standardized Dissimilarity Score (SDS)",side=2,cex=2,col="black",line=3.5)
	
	#grid:
	clip(1-0.1,length(Error_rate)+0.1,0-0.1,max(Zscore)+0.1)
  	abline(h=seq(0,2.5,by=0.5),col="gray",lty=3)
  	abline(v=seq(1,length(Error_rate),by=1),col="gray",lty=3)

  	points(Zscore, type="b",pch=16, col="black",cex=0.7)
  
  	for(i in 1:length(median)){
	    if(median[i]<0){
  				abline(	h=median[i],
  						lwd=1.5,
  						lty=2,
  						col=colors[i]
  					)	
		}else{
    		#for a rectangle:
    		#for a function y=mx+n
    	
    		if(median[i]>max(Zscore)){
    			print(	paste("SDS for ",
    					as.character(unlist(INFO[1,i]))," (sample number ",i,") "," is much higher than expected. SDS = ",median[i],".",sep=""))
    		}else{
    			y1=Zscore[max(which(Zscore<=median[i]))]
    			y2=Zscore[min(which(Zscore>=median[i]))]
    	
    			x1=max(which(Zscore<=median[i]))
    			x2=min(which(Zscore>=median[i]))
    	
    			m=(y2-y1)/(x2-x1)
    			n=y1-(m*x1)
    			s=(median[i]-n)/m
    
    			lines(x=c((1-0.4),s),y=c(median[i],median[i]),lwd=4,lty=2,col=colors[i])
    			lines(x=c(s,s),y=c(median[i],-0.1),lwd=4,lty=2,col=colors[i])
    		
    			# if(median[i]>=0.5){
    				# text(x=1,y=median[i]+0.03,paste(as.character(unlist(INFO[1,i]))," ,SDS = ",round(median[i],digits=2),sep=""),col=colors[i],pos=4)
    			# }
    		}
    	}	
	}
	##########################################################################################################################################################################
	#add 95%,5% Quantile:
  	polygon(	x=c(seq(1:5),rev(seq(1:5))),
  				y=c(konf[seq(1,10,by=2)],konf[rev(seq(2,10,by=2))]),
  				col="#99222233",
  				border=NA
  			)
	
 	#add inner legend:
  	legend(	c(1,max(Zscore)+0.1),
  			legend=paste("1000Genomes",population,"Reference Curve",sep=" "),
  			cex=2,
  			bty="n",
  			lty=1
  		)
  	#legend for samples:
  	  legend(c(1,max(Zscore)+0.1),legend=paste(rep(" ",dim(INFO)[2]),rep(" ",dim(INFO)[2]),"(",1:dim(INFO)[2],") ",as.character(unlist(INFO[1,])),sep=""),ncol=1,bty="n",xpd=T,cex=2,lty=2,col=colors,adj=0,x.intersp=0.5,pt.cex=4,lwd=4,title="",title.adj=0)


