#highlights sites-of-interest across all nucleotide sites in gene


#' Highlight sites-of-interest obtained from call to \code{\link{summary.mutPPAs}} if required.
#'   
#' @param x a \code{"summary.mutPPAs"} object.
#' @param prior a scalar used to select which results to plot according to the
#' prior probability of association. If \code{NULL} then defaults to the
#' smallest prior PA.
#' @param criteria character denoting which criteria to select sites-of-interest on. Can
#' take the values c("less", stringent").
#' @param \dots not used.
#' @author TJ McKinley
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
#' @export plot.sitesofinterest

plot.sitesofinterest<-function(x, prior = NULL, criteria = c("less", "stringent"), ...)
{
	if(class(x)!="summary.mutPPAs") stop("'x' is not a 'summary.mutPPAs' object")
	if(is.numeric(prior)==FALSE & is.null(prior)==FALSE) stop("'prior' argument is not a numeric scalar or NULL")
	if(length(prior)>1)
	{
		print("'prior' argument has length > 1, and so only first element is used")
		prior<-prior[1]
	}
	if(length(criteria) > 1) criteria <- criteria[1]
	matchtype <- match(criteria, c("less", "stringent"))
	if(is.na(matchtype)) stop("'criteria' not recognised")
	
	sites1<-x$sitesofinterest
	nsites<-nrow(sites1)
	if(is.null(nsites))
	{
		cat("No sites-of-interest\n")
		return(1)
	}
	if(nsites == 0)
	{
		cat("No sites-of-interest\n")
		return(1)
	}
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
	colnames(sites1)[1] <- "sites"
	if(criteria == "less") sites1 <- subset(sites1, select = c("sites", "Less stringent"))
	else sites1 <- subset(sites1, select = c("sites", "Stringent"))
	if(ncol(sites1) != 2) stop("Can't select on criterion for this sample")
	colnames(sites1)[2] <- "entropy"
	sites1 <- as.data.frame(sites1)
	
	#expand to whole data frame
	codons <- rep(1:3, times = x$nnuc / 3)
	sites1$codon <- codons[sites1$sites]
	sites1$codon <- factor(sites1$codon, levels = as.character(1:3))
	sites1$value <- rep(1, nrow(sites1))
	
	temp.plot <- ggplot(sites1, aes(sites, value, fill = codon)) + coord_cartesian(xlim = c(0, x$nnuc+1), ylim = c(0, 1))
#	temp.plot <- temp.plot + geom_segment(aes(x = xmin, y = gene, xend = length, yend = gene), data = codes1, arrow = arrow(angle = 90, ends = "both", length = unit(0.1, "inches")), colour = "darkgrey")
	#temp.plot <- temp.plot + geom_linerange(aes(x = sites, ymin = y, ymax = ymax))
	temp.plot <- temp.plot + geom_bar(width = 1, stat = "identity") + xlab("Nucleotide position") + ylab("") + labs(color = "Codon")
	print(temp.plot)
}
