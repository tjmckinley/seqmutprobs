\name{summary.mutPPAs}
\alias{summary.mutPPAs}
\title{Summaries of sequence information and posterior probabilities of association
from call to \code{\link{seqtoPPAs}}}
\usage{
  \method{summary}{mutPPAs} (object, thresh = 0.5,
    digits = 2, ...)
}
\arguments{
  \item{object}{a \code{"mutPPAs"} object, usually as a
  result of a call to \code{\link{seqtoPPAs}}.}

  \item{thresh}{a numerical value between 0 and 1 such that
  all sites with PPA>thresh are returned.}

  \item{digits}{a positive integer controlling how PPAs are
  rounded in output.}

  \item{\dots}{not used.}
}
\value{
  Returns a \code{"summary.mutPPAs"} object, essentially a
  list containing various summary measures of the
  corresponding \code{"mutPPAs"} object: \itemize{
  \item{nseq}{a vector containing the number of sequences
  in each sample.} \item{nnuc}{a scalar containing the
  total number of nucleotide sites in each gene segment.}
  \item{pstar}{a scalar containing the value of p*.}
  \item{samp_names}{a character vector containing the names
  of each sample.} \item{basedist}{a matrix of containing
  base distributions for each unique site.} \item{nucind}{a
  vector of length \code{nnuc} specifying which column of
  \code{basedist} corresponds to each nucleotide site.}
  \item{rem_sites}{a numeric vector containing the
  locations of any sites removed from the analysis.}
  \item{test_sites}{a numeric vector containing the
  locations of the subset of sites tested (if missing then
  all sites tested).} \item{priorPA}{a numeric vector
  containing prior probabilities of association.}
  \item{warning_sites}{a list containing information on
  sites that are valid but which have been removed due to a
  numerical precision problem (only applicable when
  \code{estimate="full"} and \code{supp_output=TRUE}).}
  \item{supp_output}{a logical recording whether the model
  outputs are suppressed.} \item{sitesofinterest}{a matrix
  containing the subset of sites-of-interest evaluated for
  the relevant criteria at a given threshold. The threshold
  is applied to the PPAs corresponding to the smallest
  prior PA.} \item{entropy}{a matrix containing values for
  the normalised continuous Kullback-Liebler divergence
  derived in McKinley et al. (2012). Comparisons are
  between the first sample and each of the subsequent
  samples, as well as the arithmetic mean of these entropy
  measures. The normalisation factor is estimated
  numerically.} \item{hyp_output}{a character vector used
  for printing.} \item{hyp_names}{a character vector used
  for printing.} \item{hyp_id}{a vector recording which
  criteria contain sites-of-interest.} \item{thresh}{a
  scalar recording the threshold applied to the PPAs.} }
}
\description{
  \code{summary} method for class \code{"mutPPAs"}
}
\details{
  Function prints some summary statistics to the screen.
}
\examples{
##read in data from fasta files
stock <- system.file("extdata/stock.fasta",
package = "seqmutprobs")
R01093seqW2 <- system.file("extdata/R01093seqW2.fasta",
package = "seqmutprobs")
R01093seqW4 <- system.file("extdata/R01093seqW4.fasta",
package = "seqmutprobs")

ref <- system.file("extdata/reference.fasta",
package = "seqmutprobs")

##combine into ordered list of 'alignment' objects
hiv_filenames <- list(stock = stock, R01093seqW2 = R01093seqW2,
R01093seqW4 = R01093seqW4)

##screen for sites-of-interest based on extracting subset of 'top' models
##and suppressing the return of model outputs for individual sites
hiv_muts <- seqtoPPAs(hiv_filenames, ref)
hiv_muts

##print summaries to screen
summary(hiv_muts, draw = TRUE)
}
\author{
  TJ McKinley
}
\references{
  McKinley et al., PLoS Comp. Biol., 7 (3), e1002027,
  (2011). doi: 10.1371/journal.pcbi.1002027
}
\seealso{
  \code{\link{seqtoPPAs}}, \code{\link{extract_site_info}},
  \code{\link{print.mutPPAs}},
  \code{\link{print.mutPPAs.list}},
  \code{\link{print.summary.mutPPAs}}
}

