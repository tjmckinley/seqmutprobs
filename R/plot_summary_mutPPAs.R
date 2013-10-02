#plot method for "summary.mutPPAs" objects


#' Plots posterior probabilities of association against entropy for those
#' sites-of-interest obtained from call to \code{\link{summary.mutPPAs}}
#' 
#' \code{plot} method for class \code{"summary.mutPPAs"}
#' 
#' Plots posterior probabilities of association against relative entropy for
#' sites-of-interest as obtained from a call to \code{summary.mutPPAs}. If
#' plotting to on-screen devices (such as \code{X11} and \code{quartz}
#' devices), then it attempts to set an optimum plot width and height for
#' visualisation, else this must be set manually.
#' 
#' @param x a \code{"summary.mutPPAs"} object.
#' @param prior a scalar used to select which results to plot according to the
#' prior probability of association. If \code{NULL} then defaults to the
#' smallest prior PA.
#' @param entropy a character corresponding to whether to plot the "max" or
#' "mean" of the absolute relative entropy values.
#' @param \dots not used.
#' @author TJ McKinley
#' @seealso \code{\link{seqtoPPAs}}, \code{\link{print.mutPPAs}},
#' \code{\link{print.mutPPAs.list}}, \code{\link{print.summary.mutPPAs}},
#' \code{\link{summary.mutPPAs}}
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
#' 
#' ##plot distributions
#' plot(summary(hiv_muts))
#' 
#' @method plot summary.mutPPAs
#' @export plot.summary.mutPPAs

plot.summary.mutPPAs<-function(x, prior=NULL, entropy=c("max","mean"), ...)
{
	if(class(x)!="summary.mutPPAs") stop("'x' is not a 'summary.mutPPAs' object")
	if(is.numeric(prior)==FALSE & is.null(prior)==FALSE) stop("'prior' argument is not a numeric scalar or NULL")
	if(length(prior)>1)
	{
		print("'prior' argument has length > 1, and so only first element is used")
		prior<-prior[1]
	}
	if(entropy[1]!="max" & entropy[1]!="mean") stop("Wrong input for 'entropy' argument")
	
	
	sites1<-x$sitesofinterest
	nsites<-nrow(sites1)
	if(nsites==0) stop("No sites-of-interest")
	if(is.null(prior))
	{
		#search for location of smallest prior
		prior<-unique(colnames(sites1)[2:ncol(sites1)])
		prior<-sapply(strsplit(prior,"\\("),function(x) x[[2]])
		prior<-sapply(strsplit(prior,"\\)"),function(x) x[[1]])
		prior1<-min(as.numeric(prior))
		prior<-paste("(",prior1,")",sep="")
		sites1<-sites1[,c(1,which(colnames(sites1)==prior))]
		if(!is.matrix(sites1)) sites1<-matrix(sites1,nrow=1)
	}
	else
	{
		prior1<-prior
		prior<-paste("(",prior,")",sep="")
		sites1<-sites1[,c(1,which(colnames(sites1)==prior))]
		if(!is.matrix(sites1)) sites1<-matrix(sites1,nrow=1)
		if(ncol(sites1)==1) stop("'prior' argument doesn't exist in 'x$sitesofinterest'")
	}
	colnames(sites1)[2:ncol(sites1)]<-x$hyp_names[x$hyp_names!=""]
	#set up colours and point symbols
	cols<-rainbow(n=nsites)
	pchs<-rep(21:25,ceiling(nsites/5))
	pchs<-pchs[1:nsites]
	ymin<-unlist(sites1[,2:ncol(sites1)])
	ymin<-ymin[!is.na(ymin)]
	ymax<-max(ymin)
	ymin<-min(ymin)
	if(entropy[1]=="max") entcol<-ncol(x$entropy)
	else entcol<-ncol(x$entropy)-1
	xmin<-unlist(x$entropy[,entcol])
	xmin<-xmin[!is.na(xmin)]
	xmax<-max(xmin)
	xmin<-min(xmin)
	#plot PPA versus entropy
	
	#calculate outer margin necessary for plot
	ncol.legend<-ceiling(nsites/25)
	
	#open new graphics device and set parameters
	if(dev.cur()==1) dev<-0
	else
	{
		dev<-names(dev.cur())
		if(length(grep("X11",dev))>0 | length(grep("quartz",dev))>0 | length(grep("windows",dev))>0) dev<-0
		else dev<-1
	}
	if(dev==0) dev.new(width=7*(ncol(sites1)-1)+0.5*ncol.legend,height=7)
	par(mfrow=c(1,ncol(sites1)-1),oma=c(0,0,2,4*ncol.legend))
	for(j in 2:ncol(sites1))
	{
		cols1<-cols[!is.na(sites1[,j])]
		pchs1<-pchs[!is.na(sites1[,j])]
		entropy1<-x$entropy[!is.na(sites1[,j]),entcol]
		sites2<-sites1[!is.na(sites1[,j]),]
		if(!is.matrix(sites2)) sites2<-matrix(sites2,nrow=1)
		plot(entropy1,sites2[,j],bg=cols1,pch=pchs1,ylab=paste("PPA (prior PA=",prior1,")",sep=""),xlab=paste(ifelse(entropy[1]=="max","Max.","Mean"),"entropy"),main=paste(colnames(sites2)[j]),ylim=c(ymin,ymax),xlim=c(xmin,xmax))
	}
	legend(par("usr")[2]*1.05,mean(par("usr")[3:4]),legend=sites1[,1],pch=pchs,pt.bg=cols,xpd=NA,yjust=0.5,ncol=ncol.legend)
	mtext(text=paste(x$genes),side=3,font=2,cex=1.5,outer=T)
}
