#plot method for "mutPPAs.list" objects


#' Plots summaries of sequence information and posterior probabilities of
#' association from call to \code{\link{seqtoPPAs}}
#' 
#' \code{plot} method for class \code{"mutPPAs.list"}
#' 
#' Plots posterior probabilities of association against relative entropy for
#' sites-of-interest as obtained from a call to \code{summary.mutPPAs}.
#' Produces one plot for each element of the \code{"summary.mutPPAs.list"}
#' object. If plotting to on-screen devices (such as \code{X11} and
#' \code{quartz} devices), then it attempts to set an optimum plot width and
#' height for each element of the object for visualisation, else it assumes a
#' fixed height and width for each element which must be set manually.
#' 
#' @param x a \code{"mutPPAs.list"} object, usually as a result of a call to
#' \code{\link{seqtoPPAs}}.
#' @param thresh a numerical value between 0 and 1 such that all sites with
#' PPA>thresh are returned.
#' @param digits a positive integer controlling how PPAs are rounded in output.
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
#' @seealso \code{\link{seqtoPPAs}}, \code{\link{summary.mutPPAs}},
#' \code{\link{print.mutPPAs}}, \code{\link{print.summary.mutPPAs}},
#' \code{\link{summary.mutPPAs.list}}
#' @references McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi:
#' 10.1371/journal.pcbi.1002027
#'
#' @method plot mutPPAs.list
#' @export plot.mutPPAs.list

plot.mutPPAs.list<-function(x,thresh=0.5,digits=2,prior=NULL,entropy=c("max","mean"),type = c("all", "ent_PPA", "shannon_ent", "prior_ent", "temporal_ent"),...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="mutPPAs.list") stop("'x' is not a 'mutPPAs.list' object")
	sapply(x,function(x) if(class(x)!="mutPPAs") stop("Elements of 'x' are not 'mutPPAs' objects"))
	plot(summary(x,thresh=thresh,digits=digits),prior=prior,entropy=entropy,type=type)
}

