#plot method for "summary.mutPPAs.list" objects


#' Plots summaries of sequence information and posterior probabilities of
#' association from call to \code{\link{seqtoPPAs}}
#' 
#' \code{plot} method for class \code{"summary.mutPPAs.list"}
#' 
#' Plots posterior probabilities of association against relative entropy for
#' sites-of-interest as obtained from a call to \code{summary.mutPPAs}.
#' Produces one plot for each element of the \code{"summary.mutPPAs.list"}
#' object. If plotting to on-screen devices (such as \code{X11} and
#' \code{quartz} devices), then it attempts to set an optimum plot width and
#' height for each element of the object for visualisation, else it assumes a
#' fixed height and width for each element which must be set manually.
#' 
#' @param x a \code{"summary.mutPPAs.list"} object.
#' @param prior a scalar used to select which results to plot according to the
#' prior probability of association. If \code{NULL} then defaults to the
#' smallest prior PA.
#' @param entropy a character corresponding to whether to plot the "max" or
#' "mean" of the absolute relative entropy values.
#' @param \dots not used.
#' @author TJ McKinley
#' @seealso \code{\link{seqtoPPAs}}, \code{\link{extract_site_info}},
#' \code{\link{summary.mutPPAs}}, \code{\link{print.mutPPAs}},
#' \code{\link{print.summary.mutPPAs}}, \code{\link{summary.mutPPAs.list}}
#' @references McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi:
#' 10.1371/journal.pcbi.1002027
#'
#' @method plot summary.mutPPAs.list
#' @export plot.summary.mutPPAs.list

plot.summary.mutPPAs.list<-function(x,prior=NULL,entropy=c("max","mean"), ...)
{
	if(missing(x)) stop("'x' argument missing")
	if(class(x)!="summary.mutPPAs.list") stop("'x' is not a 'summary.mutPPAs.list' object")
	for(i in 1:length(x)) plot(x[[i]],prior=prior,entropy=entropy)
}

