\name{print.mutPPAs}
\alias{print.mutPPAs}
\title{Prints summaries of sequence information and posterior probabilities of
association from call to \code{\link{seqtoPPAs}}}
\usage{
  \method{print}{mutPPAs} (x, thresh = 0.5, digits = 2,
    ...)
}
\arguments{
  \item{x}{a \code{"mutPPAs"} object, usually as a result
  of a call to \code{\link{seqtoPPAs}}.}

  \item{thresh}{a numerical value between 0 and 1 such that
  all sites with PPA>thresh are returned.}

  \item{digits}{a positive integer controlling how PPAs are
  rounded in output.}

  \item{\dots}{not used.}
}
\value{
  Prints the number of sequences in each sample, the
  locations of nucleotides that have been removed from the
  analysis, and the number of sites with unique base
  distributions. Also returns the locations and PPAs of all
  sites with PPA>thresh. When different prior probabilities
  of association are specified, then the threshold is
  applied to PPAs corresponding to the smallest prior PA.
}
\description{
  \code{print} method for class \code{"mutPPAs"}
}
\details{
  Function prints some summary statistics to the screen.
  Acts as a wrapper function for \code{summary.mutPPAs} and
  \code{print.summary.mutPPAs}.
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
  \code{\link{summary.mutPPAs}},
  \code{\link{print.mutPPAs.list}},
  \code{\link{print.summary.mutPPAs}}
}

