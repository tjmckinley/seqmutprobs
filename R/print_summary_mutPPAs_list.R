#print method for "summary.mutPPAs.list" objects


#' Prints summaries of sequence information and posterior probabilities of
#' association from call to \code{\link{seqtoPPAs}}
#' 
#' \code{print} method for class \code{"summary.mutPPAs.list"}
#' 
#' Function prints some summary statistics for \code{"summary.mutPPAs.list"} to
#' the screen.
#' 
#' @param x a \code{"summary.mutPPAs.list"} object.
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
#' \code{\link{print.summary.mutPPAs}}, \code{\link{summary.mutPPAs.list}}
#' @references McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi:
#' 10.1371/journal.pcbi.1002027
#'
#' @method print summary.mutPPAs.list
#' @export print.summary.mutPPAs.list

print.summary.mutPPAs.list<-function(x, ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="summary.mutPPAs.list") stop("'x' is not a 'summary.mutPPAs.list' object")
	for(i in 1:length(x)) print(x[[i]])
}

