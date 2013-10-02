\name{seqtoPPAs}
\alias{seqtoPPAs}
\alias{seqtoPPAs.default}
\title{Produces PPAs from sequence data}
\usage{
  seqtoPPAs(filenames, ref_file,
    format = c("fasta", "clustal", "phylip", "mase", "msf", "bam", "pileup"),
    ref_format = c("fasta", "clustal", "phylip", "mase", "msf"),
    estimate = c("top", "full"),
    criteria = c("both", "stringent", "less"),
    supp_output = TRUE, priorPA = c(0.001, 0.01, 0.05),
    c = 20, samp_names = NULL, pstar = NULL, sites = NA,
    genes = NA, cov_thresh = 5, nswitch_to_supp_output = 5,
    mc.cores = 1, ...)

  \method{seqtoPPAs}{default} (filenames, ref_file,
    format = c("fasta", "clustal", "phylip", "mase", "msf", "bam", "pileup"),
    ref_format = c("fasta", "clustal", "phylip", "mase", "msf"),
    estimate = c("top", "full"),
    criteria = c("both", "stringent", "less"),
    supp_output = TRUE, priorPA = c(0.001, 0.01, 0.05),
    c = 20, samp_names = NULL, pstar = NULL, sites = NA,
    genes = NA, cov_thresh = 5, nswitch_to_supp_output = 5,
    mc.cores = 1, ...)
}
\arguments{
  \item{filenames}{a list of paths to files containing
  alignments, in correct order starting with set of
  sequences for the initial biological sample.}

  \item{ref_file}{a character containing the path to a file
  containing the reference sequence(s). For NGS data this
  file must contain the same number of genes as are listed
  in the ".bam" file for the initial biological sample, and
  in the same order.}

  \item{format}{a character specifying the format for each
  input file. Must be one of "fasta", "clustal", "phylip",
  "mase", "msf", "bam" or "pileup", and all input files
  must be in the same format.}

  \item{ref_format}{a character specifying the format for
  the reference file. Must be one of "fasta", "clustal",
  "phylip", "mase" or "msf", although if input files are in
  "bam" format then \code{ref_format} must be a "fasta"
  file.}

  \item{estimate}{a character specifying whether the full
  or top model sets are required.}

  \item{criteria}{a character specifying whether the less
  stringent or stringent criteria (see \code{"Details"}),
  or both, are to be used for screening.}

  \item{supp_output}{a logical denoting whether to suppress
  the return of the individual model outputs.}

  \item{priorPA}{a numeric vector containing prior
  probabilities of association.}

  \item{c}{a scalar quantity specifying how the `top' model
  sets are selected (see \code{"Details"}).}

  \item{samp_names}{a character vector containing the names
  of each of the samples in \code{seqdata}. If \code{NULL}
  then it is generated from the samples in \code{seqdata}}

  \item{pstar}{a scalar quantity defining the
  per-nucleotide mutation probability under the random
  mutation hypothesis. If \code{NULL} then it is calculated
  directly from the samples in \code{seqdata}.}

  \item{sites}{an integer vector specifying which subset of
  nucleotide sites are to be examined. If \code{NA} then
  all sites (other than insertions and/or sites with
  invalid coverage) are evaluated.}

  \item{genes}{a vector of characters to be used with NGS
  data to specify the genes we wish to extract. Can be
  missing, in which case we compare all segments present in
  the initial sample.}

  \item{cov_thresh}{a scalar denoting the minimum coverage
  required for a site to be considered for analysis (only
  used for NGS data due to hetergeneous coverage).}

  \item{nswitch_to_supp_output}{a scalar defining the
  default number of samples above which to automatically
  suppress the model outputs.}

  \item{mc.cores}{a scalar defining the number of cores
  over which to parallelise the code. If set to \code{NA},
  then this attempts to automatically ascertain how many
  cores are available, else it uses a single core.}

  \item{\dots}{not used.}
}
\value{
  A list where each element contains an object of class
  \code{mutPPAs.list} corresponding to each gene, these are
  basically lists including elements: \itemize{
  \item{nseq:}{ if the input alignments are from capillary
  sequencing, then this is a vector the same length as
  \code{seqdata} containing the number of sequences in each
  sample. If the alignments are from NGS data, then this is
  a list where each element of the list contains the
  maximum and minimum coverage counts across the gene
  segment of interest.} \item{nnuc:}{ a scalar containing
  the total number of nucleotide sites in each gene
  segment.} \item{pstar:}{ a scalar containing the value of
  p*.} \item{samp_names:}{ a character vector the same
  length as \code{seqdata} containing the names of each
  sample.} \item{basedist:}{ a matrix of containing base
  distributions for each unique site.} \item{nucind:}{ a
  vector of length \code{nnuc} specifying which column of
  \code{basedist} corresponds to each nucleotide site.}
  \item{cons:}{ a character vector of length \code{nnuc}
  containing the consensus sequence.} \item{rem_sites:}{ a
  numeric vector containing the locations of any sites
  removed from the analysis.} \item{test_sites:}{ a numeric
  vector containing the locations of the subset of sites
  tested (if missing then all sites tested).}
  \item{priorPA:}{ a numeric vector corresponding to the
  prior probabilities of association.} \item{PPAs:}{ a list
  of length 5. If \code{output="full"} then the element
  \code{models_num} is a single matrix, and \code{hyps} is
  a list of length 2, where each element corresponds to a
  screening criteria. Otherwise each element described
  below (with the exception of \code{nprior}) is a nested
  list with (number of unique sites x length of
  \code{priorPA} x number of criteria) elements. The forms
  for each of these final elements is described below:}
  \itemize{ \item{models_num:}{ the top-level elements
  contain matrix representations of all possible models for
  all unique sites. For each model matrix, if \code{nsamp}
  is the number of samples, then each row corresponds to a
  model, and the first \code{nsamp} columns correspond to
  the choice of model for each sample, and the second
  \code{nsamp} columns correspond to a set of dependency
  indicators corresponding to whether samples are assumed
  to be drawn from exactly the same underlying
  distribution(s) or not. These are converted to a prettier
  form when printed.} \item{hyps:}{ the top-level elements
  contain binary vectors corresponding to which models in
  the corresponding element of \code{models_num} are
  consistent with the null (0) and alternative (1)
  hypotheses.} \item{nprior:}{ a numeric vector of length
  3, containing the total number of possible models, and
  the number of those corresponding to the null and
  alternative hypotheses.} \item{PPAs:}{ the top-level
  elements contain numeric vectors corresponding to the
  PPAs for each model in the corresponding element of
  \code{models_num}.} \item{norms:}{ the top-level elements
  are scalars corresponding to the normalising constants
  for each set of PPAs. (These are missing when
  \code{estimate="approx"}.)} } \item{genes:}{ a vector of
  characters recording the names of the genes.}
  \item{estimate:}{ a character recording the type of
  estimate produced.} \item{supp_output:}{ a logical
  recording whether the model outputs are suppressed.}
  \item{hyp_PPAs:}{ a nested list with (length of
  \code{priorPA} x number of criteria x number of unique
  sites) elements, where each top-level element is a scalar
  corresponding to the PPA that the given site is a
  site-of-interest given the relevant screening criteria
  (denoted PPA_SI in references).} }

  \code{NULL}
}
\description{
  Takes a list of samples, where each element of the list
  contains a set of sequences, and computes a set of
  corresponding posterior probabilities of association
}
\details{
  Function takes a list of paths-to-files where each file
  corresponds to a set of alignments from a given sample.
  These must be in the correct order for comparison (i.e.
  the initial sample must be first and so on). It is
  assumed that all sequences that will be entered will have
  been aligned to the same reference file. If
  \code{estimate="full"} then the PPAs are calculated based
  on every possible model. The memory requirements grow
  very quickly as the number of samples and/or the number
  of unique sequence sets grow (the latter in particular is
  an issue for long sequences produced by NGS
  methods)---see McKinley et al. (2012) for more details.
  If \code{estimate="top"}, then the PPAs are evaluated
  based only on a subset of the `top' models. This uses an
  efficient search algorithm defined in McKinley et al.
  (2012) to select for all models with PPAs where the ratio
  of the candidate model to the `best' model is less than
  \code{c}. If \code{supp_output=TRUE}, then none of the
  individual model outputs for each unique site are
  returned, only the overall PPAs relating to whether sites
  are sites-of-interest or not.

  The function can also read in pile-up tables produced
  using \code{samtools}. In Unix-type systems these can be
  generated using the command: \code{samtools mpileup -f
  REFERENCE BAMFILE > PILEUP.txt}, where REFERENCE is
  replaced by the path to the reference file, BAMFILE is
  replaced with the path to the corresponding BAM file, and
  PILEUP is replaced by the name of the output file. The
  reference file must be indexed in advanced (e.g. if the
  refernce file is in fasta format, then \code{samtools
  faidx ref.fasta} should do the trick. If BAM files are
  used as the input to \code{seqtoPPAs}, the function does
  these steps automatically.

  If the number of samples being evaluated is greater than
  or equal to the value in \code{nswitch_to_supp_output},
  then the function automatically suppresses the return of
  all of the individual model outputs. This is to protect
  against excess memory usage for certain situations where
  even the \code{"top"} routine can produce very large
  large numbers of models (i.e. for sites where there are
  little or no mutations).

  Due to the heterogenous coverage present in NGS data, the
  \code{cov_thresh} argument defines a minimum number of
  reads per site that must be present in the intial
  sequence in order for that site to be considered for
  analysis. Due to the large size of many BAM files, the
  \code{genes} argument allows subsets of the data to be
  analysed via gene segment.

  The screening criteria are defined in McKinley et al.
  (2010, 2012). Different prior probabilities of
  association can be specified in the \code{priorPA}
  argument, and it is also possible to screen for a subset
  of pre-defined sites using the \code{sites} argument. The
  function outputs a \code{"mutPPAs"} or
  \code{"mutPPAs.list"} object which have generic
  \code{'print'} and \code{'summary'} functions (and the
  former is also subsettable).
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

##screen for sites-of-interest based on extracting subset of 'top' models
##and returning model outputs for individual sites for quicker post-processing
hiv_muts <- seqtoPPAs(hiv_filenames, ref, supp_output = FALSE)
hiv_muts

##screen for sites-of-interest based on evaluating all possible models
##and suppressing the return of model outputs for individual sites
hiv_muts <- seqtoPPAs(hiv_filenames, ref, estimate = "full")
hiv_muts

##screen for sites-of-interest based on evaluating all possible models
##and returning model outputs for individual sites for quicker post-processing
hiv_muts <- seqtoPPAs(hiv_filenames, ref, estimate = "full", supp_output = FALSE)
hiv_muts

##subset ouput by a given site
hiv_muts[433]
}
\author{
  TJ McKinley
}
\references{
  McKinley et al., PLoS Comp. Biol., 7 (3), e1002027,
  (2011). doi: 10.1371/journal.pcbi.1002027 McKinley et al.
  (2012), in preparation.
}
\seealso{
  \code{\link{extract_site_info}},
  \code{\link{summary.mutPPAs}},
  \code{\link{print.mutPPAs}},
  \code{\link{print.mutPPAs.list}},
  \code{\link{print.summary.mutPPAs}}
}

