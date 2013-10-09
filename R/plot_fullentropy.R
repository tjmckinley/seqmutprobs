#plots entropy across all nucleotide sites in gene


#' Plots entropy across all nucleotide sites, highlighting those
#' sites-of-interest obtained from call to \code{\link{summary.mutPPAs}} if required.
#'   
#' @param x a \code{"summary.mutPPAs"} object.
#' @param highlight_sites a logical specifying whether to colour according to 
#' whether sites are classified as a site-of-interest.
#' @param prior a scalar used to select which results to plot according to the
#' prior probability of association. If \code{NULL} then defaults to the
#' smallest prior PA.
#' @param entropy a character corresponding to whether to plot the "max" or
#' "mean" of the absolute relative entropy values.
#' @param entropy_type a character denoting the type of entropy to plot. Can
#' take the values c("kl", "shannon", "klprior").
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
#' @export plot.fullentropy

plot.fullentropy <- function(x, highlight_sites = TRUE, prior = NULL, entropy = c("max", "mean"), entropy_type = c("kl", "shannon", "klprior"), criteria = c("less", "stringent"), ...)
{
	if(class(x)!="summary.mutPPAs") stop("'x' is not a 'summary.mutPPAs' object")
	if(!is.logical(highlight_sites) | length(highlight_sites) > 1) stop("'highlight_sites' is not a logical value or it's a vector")
	if(is.numeric(prior)==FALSE & is.null(prior)==FALSE) stop("'prior' argument is not a numeric scalar or NULL")
	if(length(prior)>1)
	{
		print("'prior' argument has length > 1, and so only first element is used")
		prior<-prior[1]
	}
	if(length(entropy) > 1) entropy <- entropy[1]
	if(entropy!="max" & entropy!="mean") stop("Wrong input for 'entropy' argument")
	if(!is.character(entropy_type)) stop("'entropy_type' argument is not a character vector")
	if(length(entropy_type) > 1) entropy_type <- entropy_type[1]
	matchtype <- match(entropy_type, c("kl", "shannon", "klprior"))
	if(is.na(matchtype)) stop("'entropy_type' not recognised")
	
	if(length(criteria) > 1) criteria <- criteria[1]
	matchtype <- match(criteria, c("less", "stringent"))
	if(is.na(matchtype)) stop("'criteria' not recognised")
	
	if(highlight_sites == TRUE)
	{
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
	}
	
	#generate entropies across all nucleotide sites
	entropies <- apply(x$basedist, 2, function(x, ent_type, entropy)
	{
		x <- matrix(x, nrow = 4)
		if(ncol(x) <= 1 & ent_type == "kl") 
		{
			cat("Cannot plot K-L entropy since there are < 2 samples\n")
			return(1)
		}
		if(ent_type == "kl")
		{
			ans <- numeric(ncol(x) - 1)
			#calculate entropy
			for(j in 2:ncol(x))
			{
				#produce normalised entropy
				s2 <- sum(x[, j])
				s2 <- s2 * diag(4)
				ans1 <- apply(s2, 2, function(x, dists) entropy.fn(dists[, 1], x), dists = x)
				ans[j-1] <- entropy.fn(x[, 1], x[, j]) / max(ans1)
			}
		}
		else
		{
			if(ent_type == "shannon")
			{
				ans <- apply(x, 2, shannon.fn)
				ans <- ans / (-4 * 0.25 * log(0.25))
				ans
			}
			else
			{
				ans <- apply(x, 2, function(dist)
				{
					s <- sum(dist)
					s <- s * diag(4)
					ans1 <- apply(s, 2, klprior.fn)
					ans <- klprior.fn(dist)/max(ans1)
					ans <- 1 - ans
					ans
				})
			}
		}
		if(is.null(nrow(ans))) ans <- matrix(ans, nrow = 1)
		if(entropy == "max") ans <- apply(ans, 1, max)
		else ans <- apply(ans, 1, mean)
		ans
	}, ent_type = entropy_type, entropy = entropy)
	
	#now expand to whole gene segment
	entropies <- data.frame(entropies = entropies[x$nucind], site = 1:x$nnuc)
	
	if(highlight_sites == TRUE)
	{
		entropies$group <- rep("Unselected sites", nrow(entropies))
		entropies$group[sites1$sites] <- "Site-of-interest"
		entropies$group <- factor(entropies$group, levels = c("Site-of-interest", "Unselected sites"))
		#plot entropies
		temp.plot <- ggplot(entropies, aes(site, entropies, colour = group)) + geom_point() + xlab("Nucleotide position") + ylab("Entropy")
	}
	else
	{
		#plot entropies
		temp.plot <- ggplot(entropies, aes(site, entropies, colour = entropies)) + geom_point() + xlab("Nucleotide position") + ylab("Entropy")
	}
	print(temp.plot)
	print("NEED TO CHECK K-L PRIOR IS CORRECT WAY ROUND")
}
