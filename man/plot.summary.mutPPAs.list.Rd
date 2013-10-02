\name{plot.summary.mutPPAs.list}
\alias{plot.summary.mutPPAs.list}
\title{Plots summaries of sequence information and posterior probabilities of
association from call to \code{\link{seqtoPPAs}}}
\usage{
  \method{plot}{summary.mutPPAs.list} (x, prior = NULL,
    entropy = c("max", "mean"), ...)
}
\arguments{
  \item{x}{a \code{"summary.mutPPAs.list"} object.}

  \item{prior}{a scalar used to select which results to
  plot according to the prior probability of association.
  If \code{NULL} then defaults to the smallest prior PA.}

  \item{entropy}{a character corresponding to whether to
  plot the "max" or "mean" of the absolute relative entropy
  values.}

  \item{\dots}{not used.}
}
\description{
  \code{plot} method for class
  \code{"summary.mutPPAs.list"}
}
\details{
  Plots posterior probabilities of association against
  relative entropy for sites-of-interest as obtained from a
  call to \code{summary.mutPPAs}. Produces one plot for
  each element of the \code{"summary.mutPPAs.list"} object.
  If plotting to on-screen devices (such as \code{X11} and
  \code{quartz} devices), then it attempts to set an
  optimum plot width and height for each element of the
  object for visualisation, else it assumes a fixed height
  and width for each element which must be set manually.
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
  \code{\link{print.mutPPAs}},
  \code{\link{print.summary.mutPPAs}},
  \code{\link{summary.mutPPAs.list}}
}

