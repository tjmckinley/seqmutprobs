#generate consensus sequence from a set of samples
calc_cons<-function(seq_matrix,reference)
{
	cons<-sapply(as.list(1:ncol(seq_matrix)),function(i,x,ref)
	{
		x<-x[,i]
		y<-table(x)
		y<-y[!is.na(y)]
		if(length(y)!=0)
		{
			x<-names(y)[y==max(y)]
			if(length(x)>1)
			{
				if(is.na(match(ref[i],x)))
				{
					cat(paste("Can't generate a consensus for site,",i,"\nand can't match to reference in a sensible manner,\nso this site has been removed:\nREF = ",ref[i],"\n"))
					y<-matrix(y[1:4],ncol=1)
					rownames(y)<-c("a","c","g","t")
					cat(print(y))
					cat("\n")
					x<-NA
				}
				else x<-ref[i]
			}
		}
		else x<-NA
		x
	},x=seq_matrix,ref=reference)
	cons
}

