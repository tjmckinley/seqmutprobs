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
#' @param coverage a logical indicating whether coverage is to be plotted or not.
#' @param log_coverage a logical indicating whether coverage is to be plotted on the log-scale or not.
#' @param nrow_legend a numeric indicating the maximum number of rows in the legend 
#' (automatically updating the columns accordingly).
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

plot.fullentropy <- function(x, highlight_sites = TRUE, prior = NULL, entropy = c("max", "mean"), entropy_type = c("kl", "shannon", "klprior"), criteria = c("less", "stringent"), coverage = TRUE, log_coverage = FALSE, ncol.legend = 7, ...)
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
	
	if(length(coverage) > 1) coverage <- coverage[1]
	if(!is.logical(coverage)) stop("'coverage' is not a logical value")
	
	if(length(log_coverage) > 1) log_coverage <- log_coverage[1]
	if(!is.logical(log_coverage)) stop("'log_coverage' is not a logical value")
	
	if(length(ncol.legend) > 1) ncol.legend <- ncol.legend[1]
	if(!is.numeric(ncol.legend)) stop("'ncol.legend' is not a number")
	
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
		colnames(sites1)[2] <- "PPA"
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
	xnucind <- x$nucind
	xnucind[is.na(xnucind)] <- 1
	entropies <- data.frame(entropies = entropies[xnucind], site = 1:x$nnuc)
	entropies$entropies[is.na(x$nucind)] <- NA
	
	#generate correct labelling for y-axis
	ylabel <- switch(entropy_type, kl = "K-L entropy (relative to first sample)", shannon = "Shannon entropy", klprior = "K-L entropy (relative to prior)")
	ylabel <- ifelse(entropy == "max", paste("Max.", ylabel), paste("Mean", ylabel))
	
	if(highlight_sites == TRUE)
	{ 
		#save maximum entropy value to re-scale coverage if required
		maxent <- max(na.omit(entropies$entropies))
		
		#remove sites-of-interest from original data
		sites2 <- entropies[sites1$sites, ]
		entropies <- entropies[-sites1$sites, ]
		
		#set up coverage if required
		if(coverage == TRUE)
		{
			#calculate coverage for each sample
			coverage <- apply(x$basedist, 2, function(x)
			{
				x <- matrix(x, nrow = 4)
				x <- apply(x, 2, sum)
				x
			})
			maxcov <- max(coverage)
			if(log_coverage == TRUE) coverage <- log(coverage)
			#rescale to match entropy scale
			coverage <- maxent * (coverage / ifelse(log_coverage == TRUE, log(maxcov), maxcov))
			#expand to full segment
			coverage <- coverage[, xnucind]
			coverage[, is.na(x$nucind)] <- rep(NA, nrow(coverage))
			#collapse to data frame (calling coverage column 'entropies' for ease-of-plotting)
			coverage <- data.frame(entropies = matrix(t(coverage), ncol = 1), site = rep(1:x$nnuc, times = nrow(coverage)), sample = rep(1:nrow(coverage), each = x$nnuc))
			coverage$sample <- factor(coverage$sample)
			
			#remove missing values
			entropies <- na.omit(entropies)
			
			#set up plot with coverage layered in the background
			temp.plot <- ggplot(coverage, colour = "grey") + geom_ribbon(aes(x = site, ymax = entropies, group = sample), ymin = 0, alpha = 0.25) + xlab("Nucleotide position") + ylab(ylabel)
			#add sites
			temp.plot <- temp.plot + geom_point(aes(site, entropies), data = entropies, colour = "grey")
			#add sites-of-interest
			temp.plot <- temp.plot + geom_point(aes(site, entropies, fill = factor(site)), data = sites2, pch = 21) + guides(col = guide_legend(ncol = ceiling(nrow(sites2) / ncol.legend))) + labs(colour = "Site")
			
			#add annotation mapping maximum coverage
			temp.plot <- temp.plot + geom_text(aes(x, y, label = caption), data = data.frame(x = x$nnuc * 0.9, y = maxent, caption = ifelse(log_coverage == TRUE, paste0("Max. coverage\nlog(", maxcov,") = ", round(log(maxcov), digits = 2)), paste0("Max coverage\n", maxcov))), hjust = 0.5, vjust = 1, size = 4)
		}
		else
		{
			#remove missing values
			entropies <- na.omit(entropies)
			#plot entropies
			temp.plot <- ggplot(entropies, aes(site, entropies), colour = "grey") + geom_point() + xlab("Nucleotide position") + ylab(ylabel)
			#add sites-of-interest
			temp.plot <- temp.plot + geom_point(aes(fill = factor(site)), data = sites2, pch = 21) + guides(col = guide_legend(ncol = ceiling(nrow(sites2) / ncol.legend))) + labs(colour = "Site")
		}
	}
	else
	{
		#set up coverage if required
		if(coverage == TRUE)
		{
			#save maximum entropy value to re-scale coverage if required
			maxent <- max(na.omit(entropies$entropies))
			
			#calculate coverage for each sample
			coverage <- apply(x$basedist, 2, function(x)
			{
				x <- matrix(x, nrow = 4)
				x <- apply(x, 2, sum)
				x
			})
			maxcov <- max(coverage)
			if(log_coverage == TRUE) coverage <- log(coverage)
			#rescale to match entropy scale
			coverage <- maxent * (coverage / ifelse(log_coverage == TRUE, log(maxcov), maxcov))
			#expand to full segment
			coverage <- coverage[, xnucind]
			coverage[, is.na(x$nucind)] <- rep(NA, nrow(coverage))
			#collapse to data frame (calling coverage column 'entropies' for ease-of-plotting)
			coverage <- data.frame(entropies = matrix(t(coverage), ncol = 1), site = rep(1:x$nnuc, times = nrow(coverage)), sample = rep(1:nrow(coverage), each = x$nnuc))
			coverage$sample <- factor(coverage$sample)
						
			#remove missing values
			entropies <- na.omit(entropies)
			
			#set up plot with coverage layered in the background
			temp.plot <- ggplot(coverage, colour = "grey") + geom_ribbon(aes(x = site, ymax = entropies, group = sample), ymin = 0, alpha = 0.25) + xlab("Nucleotide position") + ylab(ylabel)
			#add sites
			temp.plot <- temp.plot + geom_point(aes(site, entropies, colour = entropies), data = entropies)
			#add annotation mapping maximum coverage
			temp.plot <- temp.plot + geom_text(aes(x, y, label = caption), data = data.frame(x = x$nnuc * 0.9, y = maxent, caption = ifelse(log_coverage == TRUE, paste0("Max. coverage\nlog(", maxcov,") = ", round(log(maxcov), digits = 2)), paste0("Max coverage\n", maxcov))), hjust = 0.5, vjust = 1, size = 4)
		}
		else
		{
			#remove missing values
			entropies <- na.omit(entropies)
			#plot entropies
			temp.plot <- ggplot(entropies, aes(site, entropies, colour = entropies)) + geom_point() + xlab("Nucleotide position") + ylab(ylabel)
		}
	}
	print(temp.plot)
	print("NEED TO CHECK K-L PRIOR IS CORRECT WAY ROUND")
}
