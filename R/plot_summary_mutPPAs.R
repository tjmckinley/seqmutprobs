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
#' @param type a character vector denoting the type of plots to output. Takes any of the
#' values c("all", "ent_PPA", "shannon_ent", "prior_ent", "temporal_ent"). If any of the
#' values is "all", then the other values will be ignored.
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

plot.summary.mutPPAs<-function(x, prior=NULL, entropy=c("max","mean"), type = c("all", "ent_PPA", "shannon_ent", "prior_ent", "temporal_ent"), ...)
{
	if(class(x)!="summary.mutPPAs") stop("'x' is not a 'summary.mutPPAs' object")
	if(is.numeric(prior)==FALSE & is.null(prior)==FALSE) stop("'prior' argument is not a numeric scalar or NULL")
	if(length(prior)>1)
	{
		print("'prior' argument has length > 1, and so only first element is used")
		prior<-prior[1]
	}
	if(entropy[1]!="max" & entropy[1]!="mean") stop("Wrong input for 'entropy' argument")
	if(!is.character(type)) stop("'type' argument is not a character vector")
	matchtype <- match(type, c("all", "ent_PPA", "shannon_ent", "prior_ent", "temporal_ent"))
	if(length(matchtype[is.na(matchtype)]) > 0) stop("Wrong input for 'type' argument")
	if(length(which(type == "all")) > 0) type <- c("ent_PPA", "shannon_ent", "prior_ent", "temporal_ent")
	
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
	
	#plot PPA versus entropy
	
	if(entropy[1] == "max") entcol <- ncol(x$entropy)
    else entcol <- ncol(x$entropy) - 1
    
    #set up list of plot objects
    plots.list <- list(NULL)
    plots.ind <- 1

	#generate data frame containing information to be plotted
	for(j in 2:ncol(sites1))
	{
		entropy1 <- x$entropy[!is.na(sites1[, j]), entcol]
		sites2 <- sites1[!is.na(sites1[, j]), ]
		if(!is.matrix(sites2)) sites2 <- matrix(sites2, nrow = 1)
		if(j == 2) entropy.dat <- data.frame(sites = sites2[, 1], x = entropy1, y = sites2[, j], names = colnames(sites1)[j])
		else entropy.dat <- rbind(entropy.dat, data.frame(sites = sites2[, 1], x = entropy1, y = sites2[, j], names = colnames(sites1)[j]))
	}
	entropy.dat$sites <- factor(entropy.dat$sites)
	ent.levels <- levels(entropy.dat$sites)
	ent.sites <- sites1[, 1]
	
	if(length(which(type == "ent_PPA")) > 0)
	{
		plots.list[[plots.ind]] <- qplot(x, y, data = entropy.dat, geom = "point", colour = sites, xlab = "K-L entropy (relative to first sample)", ylab = paste("PPA (prior PA=",prior1,")",sep="")) + facet_grid(names ~ .) + guides(col = guide_legend(ncol = ceiling(length(ent.levels)/7), title = "Sites"))
		plots.ind <- plots.ind + 1
	}
	
	#plot entropy vs Shannon entropy
	entropy.dat <- data.frame(sites = ent.sites, x = x$shannon[, entcol + 1], y = x$entropy[, entcol])
	entropy.dat$sites <- factor(entropy.dat$sites, levels = ent.levels)
	
	if(length(which(type == "shannon_ent")) > 0)
	{
		plots.list[[plots.ind]] <- qplot(x, y, data = entropy.dat, geom = "point", colour = sites, ylab = "K-L entropy (relative to first sample)", xlab = "Shannon entropy") + guides(col = guide_legend(ncol = ceiling(length(ent.levels)/7), title = "Sites"))
		plots.ind <- plots.ind + 1
	}
	
	#plot entropy vs Shannon entropy
	entropy.dat <- data.frame(sites = ent.sites, x = x$klprior[, entcol + 1], y = x$entropy[, entcol])
	entropy.dat$sites <- factor(entropy.dat$sites, levels = ent.levels)
	
	if(length(which(type == "prior_ent")) > 0)
	{
		plots.list[[plots.ind]] <- qplot(x, y, data = entropy.dat, geom = "point", colour = sites, ylab = "K-L entropy (relative to first sample)", xlab = "1 - K-L entropy (relative to prior)") + guides(col = guide_legend(ncol = ceiling(length(ent.levels)/7), title = "Sites"))
		plots.ind <- plots.ind + 1
	}
	
	#plot entropy over time
	if(length(which(type == "temporal_ent")) > 0)
	{
		entropy1 <- as.data.frame(x$entropy)
		entropy1 <- entropy1[, 1:(ncol(entropy1) - 2), drop = FALSE]
		entropy1 <- data.frame(entropy = unlist(entropy1), comparison = rep(colnames(entropy1), each = nrow(entropy1)))
		entropy1 <- cbind(sites = rep(sites1[, 1], times = nrow(entropy1)/nrow(sites1)), entropy1)
		entropy1$sites <- factor(as.character(entropy1$sites), levels = ent.levels)
		if(nrow(entropy1) > 2)
		{
			plots.list[[plots.ind]] <- ggplot(entropy1, aes(comparison, entropy, group = sites, colour = sites)) + geom_line() + guides(col = guide_legend(ncol = ceiling(length(ent.levels)/7), title = "Sites")) + ylab("K-L entropy (relative to first sample)") + xlab("Samples (in order)")
			plots.ind <- plots.ind + 1
		}
	}
	
	multiplot(plotlist = plots.list, cols = 2)
}
