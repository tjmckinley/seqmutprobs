\name{plot.summary.mutPPAs}
\alias{plot.summary.mutPPAs}
\title{Plots posterior probabilities of association against entropy for those
sites-of-interest obtained from call to \code{\link{summary.mutPPAs}}}
\usage{
  \method{plot}{summary.mutPPAs} (x, prior = NULL,
    entropy = c("max", "mean"),
    type = c("all", "ent_PPA", "shannon_ent", "prior_ent", "temporal_ent"),
    ...)
}
\arguments{
  \item{x}{a \code{"summary.mutPPAs"} object.}

  \item{prior}{a scalar used to select which results to
  plot according to the prior probability of association.
  If \code{NULL} then defaults to the smallest prior PA.}

  \item{entropy}{a character corresponding to whether to
  plot the "max" or "mean" of the absolute relative entropy
  values.}

  \item{type}{a character vector denoting the type of plots
  to output. Takes any of the values c("all", "ent_PPA",
  "shannon_ent", "prior_ent", "temporal_ent"). If any of
  the values is "all", then the other values will be
  ignored.}

  \item{\dots}{not used.}
}
\description{
  \code{plot} method for class \code{"summary.mutPPAs"}
}
\details{
  Plots posterior probabilities of association against
  relative entropy for sites-of-interest as obtained from a
  call to \code{summary.mutPPAs}. If plotting to on-screen
  devices (such as \code{X11} and \code{quartz} devices),
  then it attempts to set an optimum plot width and height
  for visualisation, else this must be set manually.
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

##plot distributions
plot(summary(hiv_muts))
}
\author{
  TJ McKinley
}
\references{
  McKinley et al., PLoS Comp. Biol., 7 (3), e1002027,
  (2011). doi: 10.1371/journal.pcbi.1002027
}
\seealso{
  \code{\link{seqtoPPAs}}, \code{\link{print.mutPPAs}},
  \code{\link{print.mutPPAs.list}},
  \code{\link{print.summary.mutPPAs}},
  \code{\link{summary.mutPPAs}}
}

