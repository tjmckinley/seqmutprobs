\name{plot_fullentropy}
\alias{plot_fullentropy}
\title{Plots entropy across all nucleotide sites, highlighting those
sites-of-interest obtained from call to \code{\link{summary.mutPPAs}} if required.}
\usage{
  plot_fullentropy(x, highlight_sites = TRUE, prior = NULL,
    entropy = c("max", "mean"),
    entropy_type = c("kl", "shannon", "klprior"),
    criteria = c("less", "stringent"), coverage = TRUE,
    log_coverage = FALSE, nrow_legend = 7, thresh = 0.5,
    digits = 2, ...)
}
\arguments{
  \item{x}{a \code{"mutPPAs"}, \code{"mutPPAs.list"} or
  \code{"summary.mutPPAs"} object.}

  \item{highlight_sites}{a logical specifying whether to
  colour according to whether sites are classified as a
  site-of-interest.}

  \item{prior}{a scalar used to select which results to
  plot according to the prior probability of association.
  If \code{NULL} then defaults to the smallest prior PA.}

  \item{entropy}{a character corresponding to whether to
  plot the "max" or "mean" of the absolute relative entropy
  values.}

  \item{entropy_type}{a character denoting the type of
  entropy to plot. Can take the values c("kl", "shannon",
  "klprior").}

  \item{criteria}{character denoting which criteria to
  select sites-of-interest on. Can take the values
  c("less", stringent").}

  \item{coverage}{a logical indicating whether coverage is
  to be plotted or not.}

  \item{log_coverage}{a logical indicating whether coverage
  is to be plotted on the log-scale or not.}

  \item{nrow_legend}{a numeric indicating the maximum
  number of rows in the legend (automatically updating the
  columns accordingly).}

  \item{thresh}{a numerical value between 0 and 1 such that
  all sites with PPA > thresh are returned.}

  \item{digits}{a positive integer controlling how PPAs are
  rounded in output.}

  \item{\dots}{not used.}
}
\description{
  Plots entropy across all nucleotide sites, highlighting
  those sites-of-interest obtained from call to
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

