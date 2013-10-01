\name{extract_site_info}
\alias{extract_site_info}
\title{Extract model summaries}
\description{Extracts summary information for the base distributions and the best supported models from \code{"mutPPAs"} object.}
\usage{
extract_site_info(object, site, criteria, 
	num_mod = 5, digits = 2, c = 20, draw = FALSE, ...)
}
\arguments{
\item{object}{ an object of class \code{"mutPPAs"}.}
\item{site}{ a positive integer corresponding to nucleotide site of interest.}
\item{criteria}{ a character describing the name of the screening criterion to extract from an object of class \code{"mutPPAs"}.}
\item{num_mod}{ a positive integer corresponding to the number of models to print to the screen.}
\item{digits}{ a positive integer controlling how output is rounded.}
\item{c}{ a scalar quantity specifying how the `top' model sets are selected (see \code{"Details"}).}
\item{draw}{ a logical specifying whether a plot of the distributions is to be drawn.}
\item{...}{ not used.}
}
\details{
If the \code{"mutPPAs"} object has been generated using the options \code{output="full"} or \code{output="top"}, then this function displays results already saved in the \code{"mutPPAs"} object. If the object was produced using \code{output="none"}, then this function applies the efficient search algorithm defined in McKinley et al. (2012) to select and display for all models with PPAs where the ratio of the candidate model to the `best' model is less than \code{c} for the specified site.
}
\value{
Returns the distributions of bases across all the samples for a given site along the gene segment, and also returns summary information for the \code{num_mod} models best supported by the data.
}
\references{McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi: 10.1371/journal.pcbi.1002027}
\author{TJ McKinley}
\seealso{\code{\link{seqtoPPAs}}}
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

##extract base distributions and top 5 models at site 945, for both
##stringent and less stringent screening criteria
extract_site_info(hiv_muts, 945, "less")
extract_site_info(hiv_muts, 945, "stringent")

##plot distributions for site 945
extract_site_info(hiv_muts, 945, draw = TRUE)
}
