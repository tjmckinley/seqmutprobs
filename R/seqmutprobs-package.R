#' Screens sequence data for sites-of-interest
#' 
#' This package compares sets of genetic sequences and screens for nucleotide
#' sites-of-interest, based on the methods developed in McKinley et al., PLoS
#' Comp. Biol. (2011) and McKinley et al. (2012), in preparation.
#' 
#' \tabular{ll}{ Package: \tab seqmutprobs\cr Type: \tab Package\cr Version:
#' \tab 1.0\cr Date: \tab 2012-01-19\cr License: \tab GPL-2\cr LazyLoad: \tab
#' yes\cr } This package takes sets of ordered alignments corresponding to
#' different biological samples, and screens each nucleotide site for evidence
#' of changes in the background distribution of bases over the set of samples,
#' using the methods and screening criteria developed in McKinley et al., PLoS
#' Comp. Biol., (2011) and McKinley et al. (2012), in preparation. Posterior
#' probabilities of association (PPAs) can be generated using the
#' \code{seqtoPPAs} function. The ordered alignments can be stored in various
#' formats (e.g. fasta, mase, phylip, clustal, msf or bam). These are read in
#' using either the \code{read.alignment} function found in the \code{seqinr}
#' package, or the \code{Rsamtools} package functionality dependent on the type
#' of data being used. The PPAs are then stored in a parsimonious format, from
#' which various detailed summaries of the output can be extracted.
#' 
#' @name seqmutprobs-package
#' @aliases seqmutprobs-package seqmutprobs
#' @docType package
#' @author TJ McKinley
#' 
#' Maintainer: TJ McKinley <tjm44@@cam.ac.uk>
#' @seealso \code{seqinr}, \code{Rsamtools}
#' @references McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi:
#' 10.1371/journal.pcbi.1002027. McKinley et al. (2012), in preparation.
#' @keywords package
#' @examples
#' 
#' ##read in data from fasta files
#' stock <- system.file("extdata/stock.fasta",
#' package = "seqmutprobs")
#' R01093seqW2 <- system.file("extdata/R01093seqW2.fasta",
#' package = "seqmutprobs")
#' R01093seqW4 <- system.file("extdata/R01093seqW4.fasta",
#' package = "seqmutprobs")
#' 
#' ref <- system.file("extdata/reference.fasta",
#' package = "seqmutprobs")
#' 
#' ##combine into ordered list of 'alignment' objects
#' hiv_filenames <- list(stock = stock, R01093seqW2 = R01093seqW2, 
#' R01093seqW4 = R01093seqW4)
#' 
#' ##screen for sites-of-interest based on extracting subset of 'top' models
#' ##and suppressing the return of model outputs for individual sites
#' hiv_muts <- seqtoPPAs(hiv_filenames, ref)
#' hiv_muts
#' 
#' ##extract base distributions and top 5 models at site 945, for both
#' ##stringent and less stringent screening criteria
#' extract_site_info(hiv_muts, 945, "less")
#' extract_site_info(hiv_muts, 945, "stringent")
#' 
NULL
