\name{plot_sitesofinterest}
\alias{plot_sitesofinterest}
\title{Highlight sites-of-interest obtained from call to \code{\link{summary.mutPPAs}} if required.}
\usage{
  plot_sitesofinterest(x, prior = NULL,
    criteria = c("less", "stringent"), ...)
}
\arguments{
  \item{x}{a \code{"mutPPAs"}, \code{"mutPPAs.list"} or
  \code{"summary.mutPPAs"} object.}

  \item{prior}{a scalar used to select which results to
  plot according to the prior probability of association.
  If \code{NULL} then defaults to the smallest prior PA.}

  \item{criteria}{character denoting which criteria to
  select sites-of-interest on. Can take the values
  c("less", stringent").}

  \item{\dots}{not used.}
}
\description{
  Highlight sites-of-interest obtained from call to
  \code{\link{summary.mutPPAs}} if required.
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
}
\author{
  TJ McKinley
}
\references{
  McKinley et al., PLoS Comp. Biol., 7 (3), e1002027,
  (2011). doi: 10.1371/journal.pcbi.1002027
}

