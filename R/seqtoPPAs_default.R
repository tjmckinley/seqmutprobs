#' @return \code{NULL}
#'
#' @rdname seqtoPPAs
#' @method seqtoPPAs default
#' @S3method seqtoPPAs default

# function that will take sequences directly and output 'mutPPAs' object
seqtoPPAs.default <- function(filenames, ref_file, format = c("fasta", "clustal", 
    "phylip", "mase", "msf", "bam", "pileup"), ref_format = c("fasta", "clustal", 
    "phylip", "mase", "msf"), estimate = c("top", "full"), criteria = c("both", "stringent", 
    "less"), supp_output = TRUE, priorPA = c(0.001, 0.01, 0.05), c = 20, samp_names = NULL, 
    pstar = NULL, sites = NA, genes = NA, cov_thresh = 5, nswitch_to_supp_output = 5, 
    mc.cores = 1, ...) {
    # check that inputs are in correct format
    if (missing(filenames)) 
        stop("'filenames' argument missing")
    if (!is.list(filenames)) 
        stop("'filenames' not in list format")
    # check input format
    if (format[1] != "mase" & format[1] != "clustal" & format[1] != "phylip" & format[1] != 
        "fasta" & format[1] != "msf" & format[1] != "bam" & format[1] != "pileup") 
        stop("File for initial sample not one of either 'fasta', 'mase', 'clustal', 'phylip', 'msf', 'bam' or 'pileup' formats")
    format <- format[1]
    sapply(filenames, function(x, format1) {
        if (length(x) != 1) 
            stop("Elements of 'filenames' are too long")
        if (!is.character(x)) 
            stop("Elements of 'filenames' are not characters")
        if (!file.exists(x)) 
            stop(paste(x, "does not exist!"))
        # create index files for .bam files if necessary
        if (format1 == "bam") 
            if (!file.exists(paste(x, ".bai", sep = ""))) 
                system(paste("samtools index", x))
    }, format1 = format)
    if (length(filenames) <= 1) 
        stop("Not enough samples to conduct analysis")
    if (missing(ref_file)) 
        stop("'ref_file' argument missing")
    if (!is.character(ref_file)) 
        stop("'ref_file' argument not a character")
    if (length(ref_file) > 1) 
        stop("'ref_file' argument must be of length 1")
    if (is.null(samp_names)) {
        samp_names <- names(filenames)
        if (is.null(samp_names)) {
            samp_names <- character(length(filenames))
            for (i in 1:length(filenames)) samp_names[i] <- paste("Sample", i)
        }
    }
    if (!is.character(samp_names)) 
        stop("'samp_names' is not a character vector")
    if (length(samp_names) != length(filenames)) 
        stop("Length of 'samp_names' and 'filenames' don't match")
    
    # check format for reference file
    if (ref_format[1] != "mase" & ref_format[1] != "clustal" & ref_format[1] != "phylip" & 
        ref_format[1] != "fasta" & ref_format[1] != "msf") 
        stop("File for reference sequence not one of either 'fasta', 'mase', 'clustal', 'phylip' or 'msf' formats")
    if ((format == "bam" | format == "pileup") & ref_format[1] != "fasta") 
        stop("Reference file needs to be in 'fasta' format when reading in data from BAM files or pileup tables")
    if (!is.na(genes[1]) & !is.character(genes)) 
        stop("'genes' in incorrect format")
    if (!is.na(cov_thresh[1]) & (!is.numeric(cov_thresh[1]) | length(cov_thresh) > 
        1)) 
        stop("'cov_thresh' is not a scalar quantity")
    if (!is.null(pstar[1]) & !is.numeric(pstar[1])) 
        stop("'pstar' is not a numeric value")
    if (!is.numeric(priorPA)) 
        stop("'priorPA' not a numeric vector")
    if (!is.numeric(c) | length(c) > 1) 
        stop("Wrong input for 'c'")
    if (estimate[1] != "full" & estimate[1] != "top") 
        stop("Wrong value for 'estimate'")
    if (criteria[1] != "both" & criteria[1] != "stringent" & criteria[1] != "less") 
        stop("Wrong value for 'criteria'")
    if (!is.logical(supp_output)) 
        stop("Wrong value for 'supp_output'")
    if (sum(apply(matrix(priorPA, nrow = 1), 2, function(x) ifelse(x <= 0 | x >= 
        1, 1, 0))) > 0) 
        stop("'priorPA' incorrectly specified")
    if ((!is.na(sites[1]) & !is.numeric(sites)) | (!is.na(sites[1]) & length(sites[sites < 
        0]) > 0)) 
        stop("'sites' argument incorrect")
    if (!is.numeric(nswitch_to_supp_output) | length(nswitch_to_supp_output) > 1) 
        stop("Wrong input for 'nswitch_to_supp_output'")
    if ((!is.numeric(mc.cores) & !is.na(mc.cores)) | length(mc.cores) > 1) 
        stop("Wrong input for 'mc.cores'")
    if (is.na(mc.cores)) 
        mc.cores <- multicore:::detectCores()
    if (!is.finite(mc.cores)) 
        mc.cores <- 1
    cat(paste("No. of cores set to", mc.cores, "\n"))
    
    # check whether the number of samples is greater than or equal to
    # 'nswitch_to_supp_output'
    if (length(filenames) >= nswitch_to_supp_output & supp_output == FALSE) {
        supp_output <- TRUE
        cat("Number of samples is greater than or equal to 'nswitch_to_supp_output',\nand so individual model outputs are suppressed\n")
    }
    
    # check that the reference file and initial sample matches up
    if (format != "bam" & format != "pileup") {
        ini <- ncol(as.matrix(read.alignment(filenames[[1]], format = format)))
        ini_ref <- read.alignment(ref_file, format = ref_format[1])
        if (length(ini_ref$seq) != 1) 
            stop("Reference file contains more than one sequence")
        ini_ref <- strsplit(ini_ref$seq[[1]], "")[[1]]
        ini_ref <- ini_ref[ini_ref != "\r"]
        if (ini != length(ini_ref)) 
            stop("The length of the sequences in the initial sample and the reference file do not match")
        rm(ini)
    } else {
        if (format == "bam") 
            ini <- scanBamHeader(filenames[[1]])[[1]]$targets else {
            pathtoperl <- paste(system.file(package = "seqmutprobs"), "/Perl/pileup2csv.pl", 
                sep = "")
            seqs <- mclapply(filenames, function(files) read.table(pipe(paste("perl", 
                pathtoperl, files, sep = " ")), header = TRUE, quote = "\"", sep = "\t"), 
                mc.cores = mc.cores)
            ini <- unique(as.character(seqs[[1]]$NAME))
            names(ini) <- ini
        }
        ini_ref <- read.alignment(ref_file, format = ref_format[1])
        ini_ref_nam <- ini_ref$nam
        # remove excess '\r's if necessary
        ini_ref_nam <- sapply(as.list(ini_ref_nam), function(x) paste(strsplit(x, 
            "\r")[[1]], collapse = ""))
        if (length(ini_ref$seq) != length(ini)) 
            stop("Reference file does not contain the correct number of Gene segments compared to the initial sample")
        ini_ref <- lapply(ini_ref$seq, function(x) {
            x <- strsplit(x, "")[[1]]
            x <- x[x != "\r"]
            x
        })
        # check Gene segments are in initial sample
        names(ini_ref) <- ini_ref_nam
        if (length(which(is.na(match(ini_ref_nam, names(ini))))) != 0) 
            stop("Some Gene segments in initial sample are not in reference file")
        # for(i in 1:length(ini)) if(ini_ref$nam[i]!=names(ini)[i]) stop('Names and/or
        # order of Gene segments in reference file do not match those in the initial
        # sample')
        if (format == "bam") 
            for (i in length(ini)) if (length(ini_ref[[i]]) != ini[i]) 
                stop("The length of the sequences in the initial sample and the reference file do not match")
        rm(ini, ini_ref_nam)
    }
    
    # create pile-up tables dependent on the format of the data
    seqdata <- create.pileup(filenames, format = format, pstar = pstar, genes = genes, 
        sites = sites, samp_names = samp_names, cov_thresh = cov_thresh, reference = ini_ref, 
        ref_file = ref_file, mc.cores = mc.cores)
    if (format != "bam" & format != "pileup") 
        seqdata1 <- list(seqdata) else {
        # check for invalid Gene segments
        seqdata1 <- names(seqdata[sapply(seqdata, function(x) ifelse(is.null(names(x)[1]), 
            1, 0)) == 1])
        seqdata <- seqdata[sapply(seqdata, function(x) ifelse(is.null(names(x)[1]), 
            1, 0)) == 0]
        if (length(seqdata) == 0) 
            stop("None of the Gene segments contain any valid sites") else {
            if (length(seqdata1) > 0) {
                if (length(seqdata1) > 2) {
                  cat(paste("Gene segments ", paste(paste(seqdata1[1:(length(seqdata1) - 
                    1)], collapse = ", "), seqdata1[length(seqdata1)], collapse = " and "), 
                    " do not contain any valid sites\nand have been removed.\n", 
                    sep = ""))
                } else {
                  if (length(seqdata1) == 2) 
                    cat(paste("Gene segments ", paste(seqdata1, collapse = " and "), 
                      " do not contain any valid sites\nand have been removed.\n", 
                      sep = "")) else cat(paste("Gene segment ", seqdata1, " does not contain any valid sites\nand has been removed.\n", 
                    sep = ""))
                }
            }
        }
        seqdata1 <- seqdata
    }
    
    PPAs_output <- list(length(seqdata))
    for (q in 1:length(seqdata1)) {
        # extract correct element of seqdata
        seqdata <- seqdata1[[q]]
        cat(paste("Analysing segment", names(seqdata1)[q], "\n"))
        # check estimate type to decide subsequent analyses
        if (estimate[1] == "full") {
            # check output type to decide subsequent analyses
            if (supp_output == FALSE) {
                # calculate log[P'(D|M)]s for all models
                lPDMs <- mutlPDMs_full(seqdata, mc.cores = mc.cores)
                # calculate PPAs based on lPDMs
                PPAs <- mutPPAs_full(lPDMs, priorPA, criteria, mc.cores = mc.cores)
                rm(lPDMs)
                # calculate final PPAs based on criteria
                PPAs$hyp_PPAs <- lapply(PPAs$PPAs$PPAs, function(ppas, hyps) {
                  # cycle through priors
                  hyp_ppas <- lapply(ppas, function(ppas, hyps) {
                    # cycle through criteria
                    hyp_ppas <- lapply(as.list(1:length(ppas)), function(i, ppas, 
                      hyps) {
                      hyps <- hyps[[i]]
                      ppas <- ppas[[i]]
                      hyp_ppas <- sum(ppas[hyps == 1])
                      hyp_ppas
                    }, ppas = ppas, hyps = hyps)
                    names(hyp_ppas) <- c("less", "stringent")
                    hyp_ppas
                  }, hyps = hyps)
                  hyp_ppas
                }, hyps = PPAs$PPAs$hyps)
                # convert to correct format cycle across priors
                PPAs$hyp_PPAs <- lapply(as.list(1:length(priorPA)), function(i, hyps) {
                  # cycle across criteria
                  hyps <- lapply(as.list(1:2), function(j, hyps, i) {
                    hyps <- sapply(hyps, function(hyps, i, j) hyps[[i]][[j]], i = i, 
                      j = j)
                    hyps
                  }, hyps = hyps, i = i)
                  names(hyps) <- c("less", "stringent")
                  hyps
                }, hyps = PPAs$hyp_PPAs)
                names(PPAs$hyp_PPAs) <- paste("priorPA_", priorPA, sep = "")
                PPAs$hyp_PPAs <- lapply(PPAs$hyp_PPAs, function(hyps) lapply(hyps, 
                  function(hyps) {
                    names(hyps) <- as.character(1:length(hyps))
                    hyps
                  }))
            } else {
                # manipulate the data and generate prior information
                PPAs <- mutPPAs_top(seqdata, c, priorPA, criteria, justsetup = TRUE, 
                  mc.cores = mc.cores)
                # produce normalising constant by generating (but not recording) each potential
                # model sequentially
                norms <- gen_exact_norm(PPAs$seqdata, PPAs$pstar, PPAs$lpriors, mc.cores = mc.cores)
                if (is.null(ncol(norms))) 
                  norms <- matrix(norms, ncol = 1)
                
                # extract precision indicators
                precs <- which((1:nrow(norms))%%5 == 0)
                precs <- norms[precs, ]
                if (is.null(ncol(precs))) 
                  precs <- matrix(precs, ncol = 1)
                precs <- apply(precs, 1, function(x) list(which(x == 1)))
                precs <- lapply(precs, function(x) x[[1]])
                # remove sites that require multiple precision
                # precs1<-lapply(as.list(1:length(precs)),function(i,x,prior,ppa) { x<-x[[i]]
                # prior<-prior[i] if(length(x)>0) { x<-match(ppa$nucind,x) x<-which(!is.na(x))
                # cat(paste('WARNING: for prior PA = ',prior,', sites ',paste(x,collapse=',
                # '),'\nmight have a precision issue and have NOT been evaluated\n',sep='')) }
                # else x<-NA x<-list(x,prior) },x=precs,prior=priorPA,ppa=PPAs)
                
                # now extract sites that need rerunning and rerun using full model before
                # discarding irrelevant output (INEFFICIENT BUT AVOIDS USING GMP/MPFR IN C)
                if (sum(sapply(precs, length)) != 0) {
                  # extract unique sites that require rerunning
                  precs <- do.call("c", precs)
                  precs <- unique(precs)
                  precs <- sort(precs)
                  # create a list containing these sites in order to run full model
                  seqdata_temp <- seqdata
                  seqdata_temp$seqdata <- seqdata_temp$seqdata[, precs]
                  if (length(precs) == 1) 
                    seqdata_temp$seqdata <- matrix(seqdata_temp$seqdata, ncol = 1)
                  colnames(seqdata_temp$seqdata) <- 1:ncol(seqdata_temp$seqdata)
                  seqdata_temp$rem_sites <- (1:seqdata_temp$nnuc)[is.na(match(seqdata_temp$seqind, 
                    precs))]
                  seqdata_temp$rem_sites <- sort(seqdata_temp$rem_sites)
                  seqdata_temp$seqind <- match(seqdata_temp$seqind, precs)
                  
                  # now generate full models for subset of sites
                  lPDMs <- mutlPDMs_full(seqdata_temp, mc.cores = mc.cores)
                  # calculate PPAs based on lPDMs
                  PPAs_temp <- mutPPAs_full(lPDMs, priorPA, criteria, mc.cores = mc.cores)
                  rm(lPDMs)
                  # #calculate final PPAs based on criteria
                  PPAs_temp$hyp_PPAs <- lapply(PPAs_temp$PPAs$PPAs, function(ppas, 
                    hyps) {
                    # cycle through priors
                    hyp_ppas <- lapply(ppas, function(ppas, hyps) {
                      # cycle through criteria
                      hyp_ppas <- lapply(as.list(1:length(ppas)), function(i, ppas, 
                        hyps) {
                        hyps <- hyps[[i]]
                        ppas <- ppas[[i]]
                        hyp_ppas <- sum(ppas[hyps == 1])
                        hyp_ppas
                      }, ppas = ppas, hyps = hyps)
                      names(hyp_ppas) <- c("less", "stringent")
                      hyp_ppas
                    }, hyps = hyps)
                    hyp_ppas
                  }, hyps = PPAs_temp$PPAs$hyps)
                  # convert to correct format cycle across priors
                  PPAs_temp$hyp_PPAs <- lapply(as.list(1:length(priorPA)), function(i, 
                    hyps) {
                    # cycle across criteria
                    hyps <- lapply(as.list(1:2), function(j, hyps, i) {
                      hyps <- sapply(hyps, function(hyps, i, j) hyps[[i]][[j]], i = i, 
                        j = j)
                      hyps
                    }, hyps = hyps, i = i)
                    names(hyps) <- c("less", "stringent")
                    hyps
                  }, hyps = PPAs_temp$hyp_PPAs)
                  # names(PPAs_temp$hyp_PPAs)<-paste('priorPA_',priorPA,sep='')
                  # PPAs_temp$hyp_PPAs<-lapply(PPAs$hyp_PPAs,function(hyps)
                  # lapply(hyps,function(hyps){names(hyps)<-as.character(1:length(hyps));hyps}))
                  PPAs_temp$hyp_PPAs <- do.call("rbind", lapply(PPAs_temp$hyp_PPAs, 
                    function(x) do.call("rbind", x)))
                  PPAs_temp <- PPAs_temp$hyp_PPAs
                  PPAs_temp <- apply(PPAs_temp, 2, function(x) {
                    x <- matrix(x, ncol = length(x)/2)
                    x <- rbind(matrix(NA, nrow(x), ncol(x)), x)
                    x <- as.numeric(x)
                    x
                  })
                  # now add in the missing PPAs
                  rowinds <- (1:nrow(norms))[-which((1:nrow(norms))%%5 == 0)]
                  norms[rowinds, precs] <- PPAs_temp
                  rm(PPAs_temp)
                }
                
                # now calculate PPAs for those sites not requiring multiple precision
                norms <- norms[-which((1:nrow(norms))%%5 == 0), ]
                if (is.null(ncol(norms))) 
                  norms <- matrix(norms, ncol = 1)
                # extract PPAs based on criteria
                hyp_PPAs <- apply(norms, 2, function(x, nprior) {
                  x <- matrix(x, ncol = nprior)
                  n <- nrow(x)
                  x <- x[-(1:(n/2)), ]
                  x
                }, nprior = length(priorPA))
                # extract norms
                norms <- apply(norms, 2, function(x, nprior) {
                  x <- matrix(x, ncol = nprior)
                  n <- nrow(x)
                  x <- x[1:(n/2), ]
                  x
                }, nprior = length(priorPA))
                if (criteria[1] == "less") 
                  norms[which((1:nrow(norms))%%2 == 1) + 1, ] <- NA
                if (criteria[1] == "stringent") 
                  norms[which((1:nrow(norms))%%2 == 1), ] <- NA
                if (criteria[1] == "less") 
                  hyp_PPAs[which((1:nrow(hyp_PPAs))%%2 == 1) + 1, ] <- NA
                if (criteria[1] == "stringent") 
                  hyp_PPAs[which((1:nrow(hyp_PPAs))%%2 == 1), ] <- NA
                
                # #adjust to remove those sites with a precision issue for(i in 1:length(precs))
                # hyp_PPAs[(i-1)*2+(1:2),precs[[i]]]<-NA for(i in 1:length(precs))
                # norms[(i-1)*2+(1:2),precs[[i]]]<-NA
                
                # convert to correct format for output cycle across priors
                hyp_PPAs <- lapply(as.list(1:length(priorPA)), function(i, hyps) {
                  # cycle across criteria
                  hyps <- lapply(as.list(1:2), function(j, hyps, i) {
                    hyps <- apply(hyps, 2, function(hyps, i, j) hyps[(i - 1) * 2 + 
                      j], i = i, j = j)
                    hyps
                  }, hyps = hyps, i = i)
                  names(hyps) <- c("less", "strict")
                  hyps
                }, hyps = hyp_PPAs)
                names(hyp_PPAs) <- paste("priorPA_", priorPA, sep = "")
                PPAs$PPAs <- list(nprior = PPAs$nprior, norms = norms)
                
                hyp_PPAs <- lapply(hyp_PPAs, function(hyps) lapply(hyps, function(hyps) {
                  names(hyps) <- as.character(1:length(hyps))
                  hyps
                }))
                
                # remove extraneous outputs and re-arrange to correct format
                PPAs["lpriors"] <- NULL
                PPAs["seqdata"] <- NULL
                PPAs["nprior"] <- NULL
                PPAs$estimate <- estimate[1]
                PPAs$supp_output <- supp_output
                PPAs$hyp_PPAs <- hyp_PPAs
                # PPAs$warning_sites<-precs1
            }
        } else {
            if (supp_output == FALSE) {
                # generate top models
                PPAs <- mutPPAs_top(seqdata, c, priorPA, criteria, mc.cores = mc.cores)
                # produce normalising constants and PPA_SIs
                norms <- mclapply(PPAs$PPAs, function(ppas) {
                  # cycle through priors
                  norms <- lapply(ppas, function(ppas) {
                    # cycle through criteria
                    norms <- lapply(ppas, function(ppas) {
                      if (!is.na(ppas[1])) {
                        ppa <- exp(ppas$lPPA)
                        if (length(which(!is.finite(log(ppa)))) > 0) {
                          norm <- ppa
                          prec <- 60
                          while (length(which(!is.finite(log(norm)))) > 0 & prec <= 
                            240) {
                            prec <- prec * 2
                            norm <- exp(mpfr(ppas$lPPA, prec))
                            if (length(which(!is.finite(log(norm)))) == 0) {
                              norm <- sum(norm)
                              ppa <- mpfr(ppas$lPPA, prec)
                              ppa <- ppa - log(norm)
                              ppa <- exp(ppa)
                              hyp <- which(ppas$hyps == 1)
                              if (length(hyp) > 0) {
                                ppa <- sum(ppa[hyp])
                                ppa <- as.numeric(ppa)
                              } else ppa <- as.numeric(0)
                            }
                          }
                          if (length(which(!is.finite(log(norm)))) > 0) 
                            stop("Precision issue with normalising constant")
                          norm <- NA
                          return(list(norm = norm, hyp_PPA = ppa))
                        } else {
                          norm <- sum(ppa)
                          hyp_PPA <- exp(ppas$lPPA - log(norm))
                          hyp_PPA <- sum(hyp_PPA[ppas$hyps == 1])
                          return(list(norm = norm, hyp_PPA = hyp_PPA))
                        }
                      } else return(NA)
                    })
                    norms
                  })
                  norms
                }, mc.cores = mc.cores)
                # extract PPAs in correct order
                hyp_PPAs <- lapply(as.list(1:length(norms[[1]])), function(i, norms) {
                  # cycle through criteria
                  hyp_PPAs <- lapply(as.list(1:length(norms[[1]][[1]])), function(j, 
                    norms, i) {
                    if (!is.na(norms[[1]][[i]][[j]][1])) {
                      hyp_PPAs <- sapply(as.list(1:length(norms)), function(k, norms, 
                        i, j) norms[[k]][[i]][[j]]$hyp_PPA, norms = norms, i = i, 
                        j = j)
                      names(hyp_PPAs) <- names(norms)
                    } else hyp_PPAs <- NA
                    hyp_PPAs
                  }, norms = norms, i = i)
                  hyp_PPAs
                }, norms = norms)
                names(hyp_PPAs) <- names(norms[[1]])
                names(hyp_PPAs[[1]]) <- names(norms[[1]][[1]])
                # extract norms in correct order
                norms <- lapply(norms, function(norms) {
                  # cycle through priors
                  norms <- lapply(norms, function(norms) {
                    # cycle through criteria
                    lapply(norms, function(norms) {
                      if (!is.na(norms[1])) 
                        return(norms$norm) else return(NA)
                    })
                  })
                  norms
                })
                # extract outputs in correct order
                models_num <- lapply(PPAs$PPAs, function(ppas) {
                  # cycle through priors
                  models <- lapply(ppas, function(ppas) {
                    # cycle through criteria
                    lapply(ppas, function(ppas) {
                      if (!is.na(ppas[1])) 
                        return(ppas$models_num) else return(NA)
                    })
                  })
                  models
                })
                hyps <- lapply(PPAs$PPAs, function(ppas) {
                  # cycle through priors
                  hyps <- lapply(ppas, function(ppas) {
                    # cycle through criteria
                    lapply(ppas, function(ppas) {
                      if (!is.na(ppas[1])) 
                        return(ppas$hyps) else return(NA)
                    })
                  })
                  hyps
                })
                ppas <- PPAs$PPAs
                ppas <- mclapply(as.list(1:length(ppas)), function(i, ppas, norms) {
                  # cycle through priors
                  ppas <- lapply(as.list(1:length(ppas[[1]])), function(j, ppas, 
                    norms, i) {
                    # cycle through criteria
                    lapply(as.list(1:length(ppas[[1]][[1]])), function(k, ppas, norms, 
                      i, j) {
                      if (!is.na(ppas[[i]][[j]][[k]][1])) {
                        if (is.na(norms[[i]][[j]][[k]])) {
                          norm <- exp(ppas[[i]][[j]][[k]]$lPPA)
                          prec <- 60
                          while (length(which(!is.finite(log(norm)))) > 0 & prec <= 
                            240) {
                            prec <- prec * 2
                            norm <- exp(mpfr(ppas[[i]][[j]][[k]]$lPPA, prec))
                            if (length(which(!is.finite(log(norm)))) == 0) {
                              norm <- sum(norm)
                              ppa <- mpfr(ppas[[i]][[j]][[k]]$lPPA, prec)
                              ppa <- ppa - log(norm)
                              ppa <- as.numeric(exp(ppa))
                            }
                          }
                          if (length(which(!is.finite(log(norm)))) > 0) 
                            stop("Precision issue with normalising constant")
                          return(ppa)
                        } else return(exp(ppas[[i]][[j]][[k]]$lPPA - log(norms[[i]][[j]][[k]])))
                      } else return(NA)
                    }, ppas = ppas, norms = norms, i = i, j = j)
                  }, ppas = ppas, norms = norms, i = i)
                  ppas
                }, ppas = ppas, norms = norms, mc.cores = mc.cores)
                ppas <- list(models_num = models_num, hyps = hyps, nprior = PPAs$nprior, 
                  PPAs = ppas, norms = norms)
                
                # remove extraneous outputs and re-arrange to correct format
                PPAs["lpriors"] <- NULL
                PPAs["seqdata"] <- NULL
                PPAs["nprior"] <- NULL
                PPAs$PPAs <- ppas
                PPAs$estimate <- estimate[1]
                PPAs$supp_output <- supp_output[1]
                PPAs$hyp_PPAs <- hyp_PPAs
            } else {
                PPAs <- mutPPAs_top(seqdata, c, priorPA, criteria, ppas = FALSE, 
                  mc.cores = mc.cores)
                # extract PPAs in correct order - cycle through priors
                hyp_PPAs <- lapply(as.list(1:length(PPAs$PPAs[[1]])), function(i, 
                  ppas) {
                  # cycle through criteria
                  ppas <- lapply(as.list(1:length(ppas[[1]][[i]])), function(j, ppas) {
                    # extract correct element
                    ppas <- sapply(ppas, function(ppas, i, j) ppas[[i]][[j]]$lPPA_final, 
                      i = i, j = j)
                    ppas
                  }, ppas = ppas)
                  ppas
                }, ppas = PPAs$PPAs)
                names(hyp_PPAs) <- names(PPAs$PPAs[[1]])
                names(hyp_PPAs[[1]]) <- names(PPAs$PPAs[[1]][[1]])
                
                # extract normalising constants in correct order - cycle through priors
                norms <- lapply(as.list(1:length(PPAs$PPAs[[1]])), function(i, ppas) {
                  # cycle through criteria
                  ppas <- lapply(as.list(1:length(ppas[[1]][[i]])), function(j, ppas) {
                    # extract correct element
                    ppas <- sapply(ppas, function(ppas, i, j) ppas[[i]][[j]]$norm, 
                      i = i, j = j)
                    ppas
                  }, ppas = ppas)
                  ppas
                }, ppas = PPAs$PPAs)
                norms <- do.call("rbind", lapply(norms, function(norms) do.call("rbind", 
                  norms)))
                colnames(norms) <- names(PPAs$PPAs)
                
                PPAs$PPAs <- list(nprior = PPAs$nprior, norms = norms)
                
                # remove extraneous outputs and re-arrange to correct format
                PPAs["lpriors"] <- NULL
                PPAs["seqdata"] <- NULL
                PPAs["nprior"] <- NULL
                PPAs$estimate <- estimate[1]
                PPAs$supp_output <- supp_output
                PPAs$hyp_PPAs <- hyp_PPAs
            }
        }
        # output class 'mutPPAs' object
        class(PPAs) <- "mutPPAs"
        PPAs_output[[q]] <- PPAs
    }
    if (length(PPAs_output) == 1) 
        PPAs_output <- PPAs_output[[1]] else {
        class(PPAs_output) <- "mutPPAs.list"
        names(PPAs_output) <- names(seqdata1)
    }
    PPAs_output
} 
