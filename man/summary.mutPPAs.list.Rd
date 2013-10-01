\name{summary.mutPPAs.list}
\alias{summary.mutPPAs.list}
\title{Produces summaries of sequence information and posterior probabilities of association from call to \code{\link{seqtoPPAs}}}
\description{\code{summary} method for class \code{"mutPPAs.list"}}
\usage{
\method{summary}{mutPPAs.list}(object, thresh=0.5, digits=2, ...)
}
\arguments{
\item{object}{ a \code{"mutPPAs.list"} object, usually as a result of a call to \code{\link{seqtoPPAs}}.}
\item{thresh}{ a numerical value between 0 and 1 such that all sites with PPA>thresh are returned.}
\item{digits}{ a positive integer controlling how PPAs are rounded in output.}
\item{\dots}{ not used.}
}
\details{
Produces some summary statistics. Acts as a wrapper function for \code{summary.mutPPAs}.
}
\value{
Produces the number of sequences in each sample, the locations of nucleotides that have been removed from the analysis, and the number of sites with unique base distributions. Also returns the locations and PPAs of all sites with PPA>thresh. When different prior probabilities of association are specified, then the threshold is applied to PPAs corresponding to the smallest prior PA. 
}
\references{McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi: 10.1371/journal.pcbi.1002027}
\author{TJ McKinley}
\seealso{\code{\link{seqtoPPAs}}, \code{\link{extract_site_info}}, \code{\link{summary.mutPPAs}}, \code{\link{print.mutPPAs}}, \code{\link{print.summary.mutPPAs}}, \code{\link{print.summary.mutPPAs.list}}}

