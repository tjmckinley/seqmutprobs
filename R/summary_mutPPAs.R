#' Summaries of sequence information and posterior probabilities of association
#' from call to \code{\link{seqtoPPAs}}
#' 
#' \code{summary} method for class \code{"mutPPAs"}
#' 
#' Function prints some summary statistics to the screen.
#' 
#' @param object a \code{"mutPPAs"} object, usually as a result of a call to
#' \code{\link{seqtoPPAs}}.
#' @param thresh a numerical value between 0 and 1 such that all sites with
#' PPA>thresh are returned.
#' @param digits a positive integer controlling how PPAs are rounded in output.
#' @param \dots not used.
#' @return Returns a \code{"summary.mutPPAs"} object, essentially a list
#' containing various summary measures of the corresponding \code{"mutPPAs"}
#' object:
#' \itemize{
#' \item{nseq}{a vector containing the number of sequences in each sample.}
#' \item{nnuc}{a scalar containing the total number of nucleotide sites in
#' each gene segment.}
#' \item{pstar}{a scalar containing the value of p*.}
#' \item{samp_names}{a character vector containing the names of each
#' sample.}
#' \item{basedist}{a matrix of containing base distributions for each
#' unique site.}
#' \item{nucind}{a vector of length \code{nnuc} specifying which column of
#' \code{basedist} corresponds to each nucleotide site.}
#' \item{rem_sites}{a numeric vector containing the locations of any sites
#' removed from the analysis.}
#' \item{test_sites}{a numeric vector containing the locations of the
#' subset of sites tested (if missing then all sites tested).}
#' \item{priorPA}{a numeric vector containing prior probabilities of
#' association.}
#' \item{warning_sites}{a list containing information on sites that are
#' valid but which have been removed due to a numerical precision problem (only
#' applicable when \code{estimate="full"} and \code{supp_output=TRUE}).}
#' \item{supp_output}{a logical recording whether the model outputs are
#' suppressed.}
#' \item{sitesofinterest}{a matrix containing the subset of
#' sites-of-interest evaluated for the relevant criteria at a given threshold.
#' The threshold is applied to the PPAs corresponding to the smallest prior PA.}
#' \item{entropy}{a matrix containing values for the normalised continuous
#' Kullback-Liebler divergence derived in McKinley et al. (2012). Comparisons
#' are between the first sample and each of the subsequent samples, as well as
#' the arithmetic mean of these entropy measures. The normalisation factor is
#' estimated numerically.}
#' \item{hyp_output}{a character vector used for printing.}
#' \item{hyp_names}{a character vector used for printing.}
#' \item{hyp_id}{a vector recording which criteria contain
#' sites-of-interest.}
#' \item{thresh}{a scalar recording the threshold applied to the PPAs.}
#' }
#' @author TJ McKinley
#' @seealso \code{\link{seqtoPPAs}}, \code{\link{extract_site_info}},
#' \code{\link{print.mutPPAs}}, \code{\link{print.mutPPAs.list}},
#' \code{\link{print.summary.mutPPAs}}
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
#' ##print summaries to screen
#' summary(hiv_muts, draw = TRUE)
#'
#' @method summary mutPPAs
#' @export summary.mutPPAs

#summary method for "mutPPAs" objects
summary.mutPPAs<-function(object,thresh=0.5,digits=2, ...)
{
	#check input
	x<-object
	if(missing(x)) stop("'object' argument missing")
	if(class(x)!="mutPPAs") stop("'object' is not a 'mutPPAs' object")
	if(!is.numeric(thresh)|length(thresh)>1|thresh>1|thresh<0) stop("'thresh' parameter not in correct format")
	if(!is.numeric(digits)|length(digits)>1|digits<0) stop("'digits' parameter not in correct format")
	
	hyps<-x$hyp_PPAs
	#extract criteria by prior
	hyps<-lapply(as.list(1:2),function(i,hyps) lapply(as.list(1:length(hyps)),function(j,hyps,i) hyps[[j]][[i]],hyps=hyps,i=i),hyps=hyps)
	hyps<-lapply(hyps,function(hyps) do.call("cbind",hyps))
	hyps<-lapply(hyps,function(hyps) hyps<-cbind(inds=as.numeric(rownames(hyps)),hyps))
	hyp_output<-NA
	hyp_names_orig<-c("Less stringent","Stringent")
	if(is.na(x$PPAs$nprior[2]))
	{
		hyp_output<-"LESS STRINGENT criterion was not evaluated"
		hyps<-hyps[2]
	}
	if(is.na(x$PPAs$nprior[3]))
	{
		hyp_output<-"STRINGENT criterion was not evaluated"
		hyps<-hyps[1]
	}

	#extract subset of sites by threshold
	sitesofinterest<-lapply(hyps,function(hyps,thresh) hyps[!is.nan(hyps[,2])&!is.na(hyps[,2])&hyps[,2]>=thresh,],thresh=thresh)
	sitesofinterest[sapply(sitesofinterest,length)==0]<-NA

	#set hypothesis ID for printing
	hyp_id<-(1:length(sitesofinterest))
	hyp_names<-NA
	#remove any criteria for which there are no sites of interest
	if(length(sitesofinterest[is.na(sitesofinterest)])>0)
	{
		hyp_id<-hyp_id[!is.na(sitesofinterest)]
		sitesofinterest<-sitesofinterest[!is.na(sitesofinterest)]
	}
	if(length(sitesofinterest)>0)
	{
		nhyp<-length(hyp_id)
		#expand to include all sites
		sitesofinterest<-lapply(sitesofinterest,function(sites,nucind)
		{
			if(is.null(ncol(sites))) sites<-matrix(sites,ncol=length(sites))
			for(i in 1:length(sites[,1]))
			{
				temp<-which(nucind==sites[i,1])
				if(i==1) output<-cbind(temp,matrix(rep(sites[i,2:ncol(sites)],length(temp)),ncol=(ncol(sites)-1),byrow=T))
				else output<-rbind(output,cbind(temp,matrix(rep(sites[i,2:ncol(sites)],length(temp)),ncol=(ncol(sites)-1),byrow=T)))
			}
			output
		},nucind=x$nucind)
		#sort into ascending order by site
		sitesofinterest<-lapply(sitesofinterest,function(x){x<-x[sort.list(x[,1]),];if(!is.matrix(x)) x<-matrix(x,1);x})
		#merge output from different hypotheses together in order to print
		allsites<-matrix(unique(do.call("c",lapply(sitesofinterest,function(x) x[,1]))),nrow=1)
		allsites<-sort(allsites)
		sitesofinterest<-lapply(sitesofinterest,function(sites,unisites)
		{
			output<-matrix(NA,length(unisites),ncol(sites)-1)
			output[unisites %in% sites[,1],]<-sites[,2:ncol(sites)]
			output
		},unisites=allsites)
		sitesofinterest<-do.call("cbind",sitesofinterest)
		colnames(sitesofinterest)<-rep(paste("(",x$priorPA,")",sep=""),nhyp)
		sitesofinterest<-cbind(sites=allsites,round(sitesofinterest,digits=digits))
		nprior<-length(x$priorPA)
		nrowsites<-nrow(sitesofinterest)
		ncolsites<-ncol(sitesofinterest)
	
		#sort and output top models
		sitesofinterest<-sitesofinterest[sort.list(sitesofinterest[,2],decreasing=T),]
		if(!is.matrix(sitesofinterest))
		{
			temp_names<-names(sitesofinterest)
			sitesofinterest<-matrix(sitesofinterest,1)
			colnames(sitesofinterest)<-temp_names
		}
		sitesofinterest[,2:ncol(sitesofinterest)]<-round(sitesofinterest[,2:ncol(sitesofinterest)],digits=digits)
		hyp_names<-character(nprior*length(hyp_id))
		j<-1
		for(i in 1:length(hyp_names))
		{
			if((i-1)%%nprior==0)
			{
				hyp_names[i]<-hyp_names_orig[hyp_id[j]]
				j<-j+1
			}
			else hyp_names[i]<-""
		}
		#calculate normalised entropy measure
		dists<-x$basedist[,x$nucind[sitesofinterest[,1]]]
		if(is.null(ncol(dists))) dists<-matrix(dists,ncol=1)
		entropy<-apply(dists,2,function(x)
		{
			dists<-matrix(x,nrow=4)
			ans<-numeric(ncol(dists)-1)
			#calculate entropy
			for(j in 2:ncol(dists))
			{
				#produce normalised entropy
				s2<-sum(dists[,j])
				s2<-s2*diag(4)
				ans1<-apply(s2,2,function(x,dists) entropy.fn(dists[,1],x),dists=dists)
				ans[j-1]<-entropy.fn(dists[,1],dists[,j])/max(ans1)
			}
			ans
		})
		if(nrow(dists)>8) entropy<-t(entropy)
		else entropy<-matrix(entropy,ncol=1)
		entropy<-cbind(entropy,apply(entropy,1,mean),apply(entropy,1,max))
		entropy<-round(entropy,digits=digits)
		colnames(entropy)<-c(paste(2:(ncol(entropy)-1),":1",sep=""),"Mean","Max")
	}
	else
	{
		sitesofinterest<-NA
		hyp_id<-NA
		entropy<-NA
	}

	#remove elements of x that are no longer necessary
	x$sitesofinterest<-sitesofinterest
	x$entropy<-entropy
	x$hyp_output<-hyp_output
	x$hyp_names<-hyp_names
	x$hyp_id<-hyp_id
	x$thresh<-thresh
	class(x)<-"summary.mutPPAs"
	x
}
