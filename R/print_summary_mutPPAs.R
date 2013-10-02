#print method for "summary.mutPPAs" objects


#' Prints summaries obtained from \code{"summary.mutPPAs"} object
#' 
#' \code{print} method for class \code{"summary.mutPPAs"}
#' 
#' Function prints some summary statistics to the screen.
#' 
#' @param x a \code{"summary.mutPPAs"} object.
#' @param \dots not used.
#' @return Prints the number of sequences in each sample, the locations of
#' nucleotides that have been removed from the analysis, and the number of
#' sites with unique base distributions. Also returns the locations and PPAs of
#' all sites with PPA>thresh. When different prior probabilities of association
#' are specified, then the threshold is applied to PPAs corresponding to the
#' smallest prior PA.
#' @author TJ McKinley
#' @seealso \code{\link{seqtoPPAs}}, \code{\link{extract_site_info}},
#' \code{\link{summary.mutPPAs}}, \code{\link{print.mutPPAs}},
#' \code{\link{print.mutPPAs.list}}
#' @references McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi:
#' 10.1371/journal.pcbi.1002027
#' @examples
#' 
#' ##read in data from fasta files
#' stock <- system.file("extdata/stock.fasta",
#' package = "seqmutprobs")
#' R01093seqW2 <- system.file("extdata/R01093seqW2.fasta",
#' package = "seqmutprobs")
#' R01093seqW4 <- system.file("extdata/R01093seqW4.fasta",
#' package = "seqmutprobs")
#' 
#' ref <- system.file("extdata/reference.fasta",
#' package = "seqmutprobs")
#' 
#' ##combine into ordered list of 'alignment' objects
#' hiv_filenames <- list(stock = stock, R01093seqW2 = R01093seqW2, 
#' R01093seqW4 = R01093seqW4)
#' 
#' ##screen for sites-of-interest based on extracting subset of 'top' models
#' ##and suppressing the return of model outputs for individual sites
#' hiv_muts <- seqtoPPAs(hiv_filenames, ref)
#' hiv_muts
#' 
#' ##plot distributions in addition to printing summaries
#' summary(hiv_muts)
#' 
#' @method print summary.mutPPAs
#' @export print.summary.mutPPAs

print.summary.mutPPAs<-function(x, ...)	
{
	#print summary to screen
	cat("\n################################################################################\n")
	if(x$estimate=="full") cat(paste("Results for ",x$genes," based on evaluating ALL models\n",sep=""))
	else cat(paste("Results for ",x$genes," based on evaluating TOP models\n",sep=""))
	cat("################################################################################\n\n")
	cat(paste(length(x$nseq),"sets of samples, containing sequences of length",x$nnuc,"nucleotides:\n"))
	cat("\n")
	if(is.list(x$nseq)) cat(paste("Initial sample (",x$samp_names[1],") ranges between ",x$nseq[[1]][1]," and ",x$nseq[[1]][2] ," sequences",sep=""))
	else cat(paste("Initial sample (",x$samp_names[1],"): ",x$nseq[1]," sequences",sep=""))
	cat("\n")
	for(i in 2:length(x$nseq))
	{
		if(is.list(x$nseq)) cat(paste("Sample ",i-1," (",x$samp_names[i],") ranges between ",x$nseq[[i]][1]," and ",x$nseq[[i]][2] ," sequences",sep=""))
		else cat(paste("Sample ",i-1," (",x$samp_names[i],"): ",x$nseq[i]," sequences",sep=""))
		cat("\n")
	}
	if(!is.list(x$nseq)) cat("\nWarning: some sites might contain less bases due to removal of incomplete data.\nCheck individual sites for details.\n")
	if(!is.na(x$rem_sites[1]))
	{
		test<-x$rem_sites
		if(length(test)>1)
		{
			test1<-test[2:length(test)]-test[1:(length(test)-1)]
			test1<-c(2,test1)
			test2<-NULL
			test3<-numeric(length(test))
			j<-0
			for(i in 1:length(test))
			{
				j<-ifelse(test1[i]==1,j+1,0)
				test3[i]<-j
			}
			test3<-which(test3==0)
			if(length(test3)==1) test2<-paste(min(test),"-",max(test),sep="")
			else
			{
				for(i in 2:length(test3))
				{
					if((test3[i]-test3[i-1])==1) test2<-c(test2,paste(test[test3[i-1]],sep=""))
					else test2<-c(test2,paste(test[test3[i-1]],"-",test[test3[i]-1],sep=""))
					if(i==length(test3))
					{
						if(test3[i]<length(test)) test2<-c(test2,paste(test[test3[i]],"-",test[length(test)],sep=""))
						if(test3[i]==length(test)) test2<-c(test2,paste(test[test3[i-1]],sep=""))
					}
				}
			}
			if(length(test2)>1)
			{
				test2<-as.list(test2)
				if(length(test2)>2) test2[1:(length(test2)-2)]<-lapply(test2[1:(length(test2)-2)],function(x) c(x,", "))
				test2[length(test2)-1][[1]]<-c(test2[length(test2)-1][[1]]," ")
				test2[[length(test2)]]<-c("and ",test2[[length(test2)]][1])
				if(length(test2)>6) test2[(1:length(test2))%%6==0]<-lapply(test2[(1:length(test2))%%6==0],function(x) c(x,"\n"))
				test2<-do.call("c",test2)
				test2<-paste(test2,collapse="")
			}
		}
		else test2<-as.character(test)
		cat(paste("\nPosition(s): ",test2," \nremoved either as possible insertion sites (i.e. initial sample consensus '-')\nor because they have too few samples and/or it is not possible to generate \na consensus base\n",sep=""))
	}
	if(is.na(x$test_sites[1])) cat("\n####### ALL valid sites tested #######\n")
	else
	{
		cat(paste("\n####### SUBSET of valid sites tested #######\n",sep=""))
		cat(paste(paste(x$test_sites,collapse=" "),"\n",sep=""))
	}
	cat(paste("\n",ncol(x$basedist),"/",x$nnuc," unique distributions of bases across the combined set of samples\n",sep=""))
	cat(sprintf("\np* = %.3g\n",x$pstar))
	
	#remove any criteria for which there are no sites of interest
	if(!is.na(x$hyp_id[1]))
	{
		hyp_names_temp<-c("less stringent","stringent")
		hyp_names_temp<-hyp_names_temp[-(x$hyp_id)]
		for(i in hyp_names_temp) cat(paste("\nNo nucleotide sites found that have PPAs > ",x$thresh," for ",i," criterion\n",sep=""))
		if(!is.na(x$sitesofinterest[1]))
		{
			cat(paste("\nAll nucleotide sites that have PPAs>",x$thresh," (for priorPA = ",x$priorPA[1],"):\n",sep=""))
			if(!is.na(x$hyp_output))
			{
				cat("\n############# ")
				cat(x$hyp_output)
				cat(" #############\n\n")
			}
			#print results to screen
			x$sitesofinterest<-rbind(colnames(x$sitesofinterest),x$sitesofinterest)
			x$sitesofinterest<-rbind(c("",x$hyp_names),x$sitesofinterest)
			x$sitesofinterest[is.na(x$sitesofinterest)]<-""
			write.table(format(x$sitesofinterest),quote=F,na="",row.names=F,col.names=F)
			#now print entropy results to the screen
			cat("\n")
			x$entropy<-rbind(colnames(x$entropy),x$entropy)
			x$entropy<-rbind(c("K-L Entropy",rep("",ncol(x$entropy)-1)),x$entropy)
			x$entropy[is.na(x$entropy)]<-""
			x$entropy<-cbind(x$sitesofinterest[,1],x$entropy)
			write.table(format(x$entropy),quote=F,na="",row.names=F,col.names=F)
			
			#now print Shannon entropy results to the screen
			cat("\n")
			x$shannon<-rbind(colnames(x$shannon),x$shannon)
			x$shannon<-rbind(c("Shannon entropy",rep("",ncol(x$shannon)-1)),x$shannon)
			x$shannon[is.na(x$shannon)]<-""
			x$shannon<-cbind(x$sitesofinterest[,1],x$shannon)
			write.table(format(x$shannon),quote=F,na="",row.names=F,col.names=F)
		}
	}
	else cat(paste("\nNo nucleotide sites found that have PPAs > ",x$thresh," for EITHER criterion\n",sep=""))
}

