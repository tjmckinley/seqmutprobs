# print method for 'mutPPAs.list' objects


#' Prints summaries of sequence information and posterior probabilities of
#' association from call to \code{\link{seqtoPPAs}}
#' 
#' \code{print} method for class \code{'mutPPAs.list'}
#' 
#' Function prints some summary statistics to the screen. Acts as a wrapper
#' function for \code{summary.mutPPAs.list} and
#' \code{print.summary.mutPPAs.list}.
#' 
#' @param x a \code{'mutPPAs.list'} object, usually as a result of a call to
#' \code{\link{seqtoPPAs}}.
#' @param thresh a numerical value between 0 and 1 such that all sites with
#' PPA>thresh are returned.
#' @param digits a positive integer controlling how PPAs are rounded in output.
#' @param \dots not used.
#' @return Prints the number of sequences in each sample, the locations of
#' nucleotides that have been removed from the analysis, and the number of
#' sites with unique base distributions. Also returns the locations and PPAs of
#' all sites with PPA>thresh. When different prior probabilities of association
#' are specified, then the threshold is applied to PPAs corresponding to the
#' smallest prior PA.
#' @author TJ McKinley
#' @seealso \code{\link{seqtoPPAs}}, \code{\link{summary.mutPPAs}},
#' \code{\link{print.mutPPAs}}, \code{\link{print.summary.mutPPAs}},
#' \code{\link{summary.mutPPAs.list}}
#' @references McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi:
#' 10.1371/journal.pcbi.1002027
#'
#' @method print mutPPAs.list
#' @export print.mutPPAs.list

print.mutPPAs.list <- function(x, thresh = 0.5, digits = 2, ...) {
    if (missing(x)) 
        stop("'x' argument missing")
    if (class(x) != "mutPPAs.list") 
        stop("'x' is not a 'mutPPAs.list' object")
    sapply(x, function(x) if (class(x) != "mutPPAs") 
        stop("Elements of 'x' are not 'mutPPAs' objects"))
    print(summary(x, thresh = thresh, digits = digits, ...))
}
 
