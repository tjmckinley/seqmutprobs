#' Produces summaries of sequence information and posterior probabilities of
#' association from call to \code{\link{seqtoPPAs}}
#' 
#' \code{summary} method for class \code{'mutPPAs.list'}
#' 
#' Produces some summary statistics. Acts as a wrapper function for
#' \code{summary.mutPPAs}.
#' 
#' @param object a \code{'mutPPAs.list'} object, usually as a result of a call
#' to \code{\link{seqtoPPAs}}.
#' @param thresh a numerical value between 0 and 1 such that all sites with
#' PPA>thresh are returned.
#' @param digits a positive integer controlling how PPAs are rounded in output.
#' @param \dots not used.
#' @return Produces the number of sequences in each sample, the locations of
#' nucleotides that have been removed from the analysis, and the number of
#' sites with unique base distributions. Also returns the locations and PPAs of
#' all sites with PPA>thresh. When different prior probabilities of association
#' are specified, then the threshold is applied to PPAs corresponding to the
#' smallest prior PA. If \code{draw = TRUE} then also specifies that plots of
#' the distributions of bases for all sites-of-interest are to be produced when
#' the summaries are printed.
#' @author TJ McKinley
#' @seealso \code{\link{seqtoPPAs}}, \code{\link{extract_site_info}},
#' \code{\link{summary.mutPPAs}}, \code{\link{print.mutPPAs}},
#' \code{\link{print.summary.mutPPAs}},
#' \code{\link{print.summary.mutPPAs.list}}
#' @references McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi:
#' 10.1371/journal.pcbi.1002027
#'
#' @method summary mutPPAs.list
#' @export summary.mutPPAs.list

summary.mutPPAs.list <- function(object, thresh = 0.5, digits = 2, ...) {
    if (missing(object)) 
        stop("'object' argument missing")
    if (class(object) != "mutPPAs.list") 
        stop("'object' is not a 'mutPPAs.list' object")
    sapply(object, function(x) if (class(x) != "mutPPAs") 
        stop("Elements of 'x' are not 'mutPPAs' objects"))
    object.out <- list(NULL)
    for (i in 1:length(object)) object.out[[i]] <- summary(object[[i]], thresh = thresh, 
        digits = digits, ...)
    class(object.out) <- "summary.mutPPAs.list"
    names(object.out) <- names(object)
    object.out
}
 
