# function to extract and print models relating to specific nucleotides to be run
# on 'mutPPAs.list' objects


#' Extract model summaries
#' 
#' Extracts summary information for the base distributions and the best
#' supported models from \code{'mutPPAs'} object.
#' 
#' If the \code{'mutPPAs'} object has been generated using the options
#' \code{output='full'} or \code{output='top'}, then this function displays
#' results already saved in the \code{'mutPPAs'} object. If the object was
#' produced using \code{output='none'}, then this function applies the
#' efficient search algorithm defined in McKinley et al. (2012) to select and
#' display for all models with PPAs where the ratio of the candidate model to
#' the `best' model is less than \code{c} for the specified site.
#' 
#' @param object an object of class \code{'mutPPAs'}.
#' @param site a positive integer corresponding to nucleotide site of interest.
#' @param criteria a character describing the name of the screening criterion
#' to extract from an object of class \code{'mutPPAs'}.
#' @param num_mod a positive integer corresponding to the number of models to
#' print to the screen.
#' @param digits a positive integer controlling how output is rounded.
#' @param c a scalar quantity specifying how the `top' model sets are selected
#' (see \code{'Details'}).
#' @param draw a logical specifying whether a plot of the distributions is to
#' be drawn.
#' @param ...  not used.
#' @return Returns the distributions of bases across all the samples for a
#' given site along the gene segment, and also returns summary information for
#' the \code{num_mod} models best supported by the data.
#' @author TJ McKinley
#' @seealso \code{\link{seqtoPPAs}}
#' @references McKinley et al., PLoS Comp. Biol., 7 (3), e1002027, (2011). doi:
#' 10.1371/journal.pcbi.1002027
#' @examples
#' 
#' ##read in data from fasta files
#' stock <- system.file('extdata/stock.fasta',
#' package = 'seqmutprobs')
#' R01093seqW2 <- system.file('extdata/R01093seqW2.fasta',
#' package = 'seqmutprobs')
#' R01093seqW4 <- system.file('extdata/R01093seqW4.fasta',
#' package = 'seqmutprobs')
#' 
#' ref <- system.file('extdata/reference.fasta',
#' package = 'seqmutprobs')
#' 
#' ##combine into ordered list of 'alignment' objects
#' hiv_filenames <- list(stock = stock, R01093seqW2 = R01093seqW2, 
#' R01093seqW4 = R01093seqW4)
#' 
#' ##screen for sites-of-interest based on extracting subset of 'top' models
#' ##and suppressing the return of model outputs for individual sites
#' hiv_muts <- seqtoPPAs(hiv_filenames, ref)
#' 
#' ##extract base distributions and top 5 models at site 945, for both
#' ##stringent and less stringent screening criteria
#' extract_site_info(hiv_muts, 945, 'less')
#' extract_site_info(hiv_muts, 945, 'stringent')
#' 
#' ##plot distributions for site 945
#' extract_site_info(hiv_muts, 945, draw = TRUE)
#' 
#' @export extract_site_info

extract_site_info <- function(object, site, criteria = c("less", "stringent"), num_mod = 5, 
    digits = 2, c = 20, draw = FALSE, ...) {
    # check that inputs are in correct format
    if (missing(object)) 
        stop("'object' argument missing")
    if (missing(site)) 
        stop("'site' argument missing")
    if (class(object) != "mutPPAs") 
        stop("'object' is not a 'mutPPAs' object")
    
    hyp_name <- criteria
    if (!is.character(hyp_name[1])) 
        stop("'hyp_name' is not a character")
    if (hyp_name[1] != "less" & hyp_name[1] != "stringent") 
        stop("'hyp_name is not either 'less' or 'stringent'")
    if (length(site) > 1) 
        stop("'site' is not a scalar")
    if (!is.numeric(site) | site < 1) 
        stop("'site' argument is incorrect")
    if (!is.numeric(num_mod) | length(num_mod) > 1) 
        stop("'num_mod' argument is wrong")
    if (!is.numeric(digits) | length(digits) > 1 | digits < 0) 
        stop("'digits' parameter not in correct format")
    if (!is.logical(draw)) 
        stop("'draw' not a logical")
    
    hyp_name_orig <- c("Less stringent", "Stringent")
    
    # for each specified nucleotide site extract and print relevant information
    for (i in site) {
        ind <- object$nucind[i]
        if (is.na(ind)) 
            cat(paste("Site", i, "was removed from analysis and/or not evaluated\n")) else {
            if (object$supp_output == FALSE) {
                if (object$estimate == "full") {
                  if (hyp_name[1] == "less") {
                    # check whether criteria is missing or not
                    if (!is.na(object$PPAs$nprior[2])) {
                      n_samp_mod <- ncol(object$PPAs$models_num)
                      output <- cbind(object$PPAs$models_num, sapply(object$PPAs$PPAs[[ind]], 
                        function(ppas, crit) ppas[[crit]], crit = 1), object$PPAs$hyps[[1]])
                      models <- output[, 1:n_samp_mod]
                      output <- output[, (n_samp_mod + 1):ncol(output)]
                      # now generate character representation of models
                      models <- apply(models, 1, function(mod, nsamp) numtochar(mod, 
                        nsamp), nsamp = ncol(models)/2)
                      # append models to output data frame
                      output <- data.frame(models, output)
                      colnames(output) <- c("models", paste("(", object$priorPA, 
                        ")", sep = ""), "Null/Alt")
                    } else stop(paste("Criteria '", hyp_name[1], "' was not evaluated in this case", 
                      sep = ""))
                  } else {
                    # check whether criteria is missing or not
                    if (!is.na(object$PPAs$nprior[3])) {
                      n_samp_mod <- ncol(object$PPAs$models_num)
                      output <- cbind(object$PPAs$models_num, sapply(object$PPAs$PPAs[[ind]], 
                        function(ppas, crit) ppas[[crit]], crit = 2), object$PPAs$hyps[[2]])
                      models <- output[, 1:n_samp_mod]
                      output <- output[, (n_samp_mod + 1):ncol(output)]
                      # now generate character representation of models
                      models <- apply(models, 1, function(mod, nsamp) numtochar(mod, 
                        nsamp), nsamp = ncol(models)/2)
                      # append models to output data frame
                      output <- data.frame(models, output)
                      colnames(output) <- c("models", paste("(", object$priorPA, 
                        ")", sep = ""), "Null/Alt")
                    } else stop(paste("Criteria '", hyp_name[1], "' was not evaluated in this case", 
                      sep = ""))
                  }
                } else {
                  if (hyp_name[1] == "less") {
                    # check whether criteria is missing or not
                    if (!is.na(object$PPAs$nprior[2])) {
                      # now extract the relevant models from the 'mutPPAs' object
                      crit <- 1
                      n_samp_mod <- ncol(object$PPAs$models_num[[ind]][[1]][[crit]])
                      models_num <- lapply(object$PPAs$models_num[[ind]], function(x, 
                        crit) x[[crit]], crit = crit)
                      hyps <- lapply(object$PPAs$hyps[[ind]], function(x, crit) x[[crit]], 
                        crit = crit)
                      ppas <- lapply(object$PPAs$PPAs[[ind]], function(x, crit) x[[crit]], 
                        crit = crit)
                      # sort into the correct order and amalgamate for printing
                      models_num <- lapply(models_num, function(models) {
                        # generate character representation of models
                        models <- apply(models, 1, function(mod, nsamp) numtochar(mod, 
                          nsamp), nsamp = ncol(models)/2)
                        models
                      })
                      if (length(models_num) > 1) {
                        # match up models
                        models <- do.call("c", models_num)
                        hyps <- do.call("c", hyps)
                        hyps <- hyps[!duplicated(models)]
                        models <- models[!duplicated(models)]
                        models_num <- lapply(models_num, function(mods, compmods) match(compmods, 
                          mods), compmods = models)
                        # now expand PPAs for each prior separately and fill in the blanks
                        ppas <- lapply(as.list(1:length(ppas)), function(i, ppas, 
                          orders) {
                          ppas <- ppas[[i]]
                          orders <- orders[[i]]
                          ppas <- ppas[orders[!is.na(orders)]]
                          orders <- match(orders, orders[!is.na(orders)])
                          ppas1 <- rep(NA, length(orders))
                          ppas1[orders[!is.na(orders)]] <- ppas
                          ppas1
                        }, ppas = ppas, orders = models_num)
                        models_num <- models
                      }
                      # append models to output data frame
                      output <- data.frame(models_num, ppas, hyps)
                      colnames(output) <- c("models", paste("(", object$priorPA, 
                        ")", sep = ""), "Null/Alt")
                    } else stop(paste("Criteria '", hyp_name[1], "' was not evaluated in this case", 
                      sep = ""))
                  } else {
                    # check whether criteria is missing or not
                    if (!is.na(object$PPAs$nprior[3])) {
                      # now extract the relevant models from the 'mutPPAs' object
                      crit <- 2
                      n_samp_mod <- ncol(object$PPAs$models_num[[ind]][[1]][[crit]])
                      models_num <- lapply(object$PPAs$models_num[[ind]], function(x, 
                        crit) x[[crit]], crit = crit)
                      hyps <- lapply(object$PPAs$hyps[[ind]], function(x, crit) x[[crit]], 
                        crit = crit)
                      ppas <- lapply(object$PPAs$PPAs[[ind]], function(x, crit) x[[crit]], 
                        crit = crit)
                      # sort into the correct order and amalgamate for printing
                      models_num <- lapply(models_num, function(models) {
                        # generate character representation of models
                        models <- apply(models, 1, function(mod, nsamp) numtochar(mod, 
                          nsamp), nsamp = ncol(models)/2)
                        models
                      })
                      if (length(models_num) > 1) {
                        # match up models
                        models <- do.call("c", models_num)
                        hyps <- do.call("c", hyps)
                        hyps <- hyps[!duplicated(models)]
                        models <- models[!duplicated(models)]
                        models_num <- lapply(models_num, function(mods, compmods) match(compmods, 
                          mods), compmods = models)
                        # now expand PPAs for each prior separately and fill in the blanks
                        ppas <- lapply(as.list(1:length(ppas)), function(i, ppas, 
                          orders) {
                          ppas <- ppas[[i]]
                          orders <- orders[[i]]
                          ppas <- ppas[orders[!is.na(orders)]]
                          orders <- match(orders, orders[!is.na(orders)])
                          ppas1 <- rep(NA, length(orders))
                          ppas1[orders[!is.na(orders)]] <- ppas
                          ppas1
                        }, ppas = ppas, orders = models_num)
                        models_num <- models
                      }
                      # append models to output data frame
                      output <- data.frame(models_num, ppas, hyps)
                      colnames(output) <- c("models", paste("(", object$priorPA, 
                        ")", sep = ""), "Null/Alt")
                    } else stop(paste("Criteria '", hyp_name[1], "' was not evaluated in this case", 
                      sep = ""))
                  }
                }
            } else {
                # produce set of 'top' models for printing
                if (hyp_name[1] == "less") {
                  # check whether criteria is missing or not
                  if (!is.na(object$PPAs$nprior[2])) {
                    # now search for 'top' models
                    seqdata <- matrix(object$basedist[, ind], ncol = 1)
                    
                    # calculate models
                    nsamp <- nrow(seqdata)/4
                    if (nsamp == 1) 
                      stop("Need more than one sample")
                    
                    priorPA <- object$priorPA
                    
                    # produce prior specifications by generating (but not recording) each potential
                    # model sequentially
                    structure <- .Call("genmodels_priors", nsamp, PACKAGE = "seqmutprobs")
                    nstructure <- length(structure)
                    totmods <- structure[nstructure - 2]
                    nalt_ls <- structure[nstructure - 1]
                    nalt_s <- structure[nstructure]
                    nnull_ls <- totmods - nalt_ls
                    nnull_s <- totmods - nalt_s
                    structure <- structure[1:(nstructure - 3)]
                    # generate priorPAs for different hypotheses
                    lpriornull_ls <- (1 - priorPA)/nnull_ls
                    lprioralt_ls <- priorPA/nalt_ls
                    lpriornull_s <- (1 - priorPA)/nnull_s
                    lprioralt_s <- priorPA/nalt_s
                    lpriors <- cbind(lpriornull_ls, lprioralt_ls, lpriornull_s, lprioralt_s)
                    lpriors <- log(lpriors)
                    
                    if (hyp_name[1] == "less") 
                      nalt_s <- NA
                    if (hyp_name[1] == "stringent") 
                      nalt_ls <- NA
                    nprior <- c(nmod = totmods, nalt_less = nalt_ls, nalt_string = nalt_s)
                    # calculate log-threshold
                    logc <- log(c)
                    
                    # generate total number of columns required to store intermediate calcs
                    ntotcol <- sum(choose(nsamp, 1:nsamp))
                    
                    PPA_mat <- apply(seqdata, 2, function(seqs, nsamp, pstar, ntotcol, 
                      logc, lpriors, structure, criteria) {
                      # generate (10 x ncol)-matrix of intermediate values for calculating PPAs
                      lPDM_int_mat <- .C("calc_lPDM_int_fn", as.integer(nsamp), as.integer(seqs), 
                        as.double(pstar), lPDM_int_mat = as.double(numeric(10 * ntotcol)), 
                        PACKAGE = "seqmutprobs")$lPDM_int_mat
                      lPPAs <- apply(lpriors, 1, function(x, seqs, nsamp, pstar, 
                        ntotcol, logc, structure, criteria, lPDM_int_mat) {
                        # sort out and remove duplicates in 'lPDM_int_mat'
                        lPDM_int_mat1 <- lPDM_int_mat
                        lPDM_int_mat <- matrix(lPDM_int_mat, 10)
                        uni <- matrix(0, 10, nsamp * 10)
                        uni_ind <- matrix(0, 10, nsamp)
                        for (i in 1:nsamp) {
                          z <- unique(lPDM_int_mat[duplicated(lPDM_int_mat[, i]), 
                            i])
                          if (length(z) > 0) {
                            for (j in z) {
                              y <- which(lPDM_int_mat[, i] == j)
                              y1 <- y[1]
                              y <- y[2:length(y)]
                              lPDM_int_mat[y, i] <- (-1e+10)
                              uni[y1, (i - 1) * 10 + (1:length(y))] <- y - 1
                              uni_ind[y1, i] <- length(y)
                            }
                          }
                        }
                        lPDM_int_mat <- as.numeric(lPDM_int_mat)
                        uni <- as.numeric(uni)
                        uni_ind <- as.numeric(uni_ind)
                        if (criteria == "both" | criteria == "less") {
                          # now calculate PPAs according to approximation routine for LESS-STRINGENT
                          # criteria
                          lPPA_mat_ls <- .Call("calc_PPAs_approx_fn", nsamp, ntotcol, 
                            logc, x[1], x[2], lPDM_int_mat1, 0, structure, length(structure)/nsamp, 
                            1, uni, uni_ind, PACKAGE = "seqmutprobs")
                          totmods <- lPPA_mat_ls[length(lPPA_mat_ls)]
                          lPPA_mat_ls <- lPPA_mat_ls[1:(length(lPPA_mat_ls) - 1)]
                          models_num <- lPPA_mat_ls[1:(2 * nsamp * totmods)]
                          models_num <- matrix(models_num, nrow = totmods, byrow = T)
                          hyp <- lPPA_mat_ls[(2 * nsamp * totmods + 1):length(lPPA_mat_ls)]
                          lPPA_mat_ls <- hyp[(totmods + 1):length(hyp)]
                          multfact <- lPPA_mat_ls[(totmods + 1):length(lPPA_mat_ls)]
                          hyp <- hyp[1:totmods]
                          lPPA_mat_ls <- lPPA_mat_ls[1:totmods]
                          lPPA_mat_ls <- list(models_num = models_num, hyps = hyp, 
                            lPPA = lPPA_mat_ls)
                        }
                        if (criteria == "both" | criteria == "stringent") {
                          # now calculate PPAs according to approximation routine for LESS-STRINGENT
                          # criteria
                          lPPA_mat_s <- .Call("calc_PPAs_approx_fn", nsamp, ntotcol, 
                            logc, x[3], x[4], lPDM_int_mat1, 1, structure, length(structure)/nsamp, 
                            1, uni, uni_ind, PACKAGE = "seqmutprobs")
                          totmods <- lPPA_mat_s[length(lPPA_mat_s)]
                          lPPA_mat_s <- lPPA_mat_s[1:(length(lPPA_mat_s) - 1)]
                          models_num <- lPPA_mat_s[1:(2 * nsamp * totmods)]
                          models_num <- matrix(models_num, nrow = totmods, byrow = T)
                          hyp <- lPPA_mat_s[(2 * nsamp * totmods + 1):length(lPPA_mat_s)]
                          lPPA_mat_s <- hyp[(totmods + 1):length(hyp)]
                          multfact <- lPPA_mat_s[(totmods + 1):length(lPPA_mat_s)]
                          hyp <- hyp[1:totmods]
                          lPPA_mat_s <- lPPA_mat_s[1:totmods]
                          lPPA_mat_s <- list(models_num = models_num, hyps = hyp, 
                            lPPA = lPPA_mat_s)
                        }
                        if (criteria == "less") 
                          lPPA_mat_s <- NA
                        if (criteria == "stringent") 
                          lPPA_mat_ls <- NA
                        # output lists
                        list(less = lPPA_mat_ls, stringent = lPPA_mat_s)
                      }, seqs = seqs, nsamp = nsamp, pstar = pstar, ntotcol = ntotcol, 
                        logc = logc, structure = structure, criteria = criteria, 
                        lPDM_int_mat = lPDM_int_mat)
                      names(lPPAs) <- paste("priorPA_", priorPA, sep = "")
                      lPPAs
                    }, nsamp = nsamp, pstar = object$pstar, ntotcol = ntotcol, logc = logc, 
                      lpriors = lpriors, structure = structure, criteria = hyp_name[1])
                    # extract correct values from 'norms' if required
                    norms <- matrix(object$PPAs$norms[, ind], ncol = 1)
                    # calculate subset of top models with exact PPAs
                    PPAs.adjust <- lapply(as.list(1:ncol(norms)), function(i, ppas, 
                      norms, names.priors, estimate, crit) {
                      # extract correct sequence
                      ppas <- ppas[[i]]
                      norms <- matrix(norms[, i], nrow = 2)
                      # cycle through priors
                      ppas <- lapply(as.list(1:ncol(norms)), function(j, ppas, norms, 
                        estimate, crit) {
                        ppas <- ppas[[j]]
                        norms <- norms[crit, j]
                        # cycle through criteria
                        if (is.na(norms)) {
                          if (estimate == "top") {
                            ppa <- exp(ppas[[crit]]$lPPA)
                            if (length(which(!is.finite(log(ppa)))) > 0) {
                              norm <- ppa
                              prec <- 60
                              while (length(which(!is.finite(log(norm)))) > 0 & prec <= 
                                240) {
                                prec <- prec * 2
                                norm <- exp(mpfr(ppas[[crit]]$lPPA, prec))
                                if (length(which(!is.finite(log(norm)))) == 0) {
                                  norm <- sum(norm)
                                  ppa <- mpfr(ppas[[crit]]$lPPA, prec)
                                  ppa <- ppa - log(norm)
                                  ppa <- as.numeric(exp(ppa))
                                }
                              }
                              if (length(which(!is.finite(log(norm)))) > 0) 
                                stop("Precision issue with normalising constant")
                              ppas[[crit]]$PPAs <- ppa
                              ppas[[crit]]$norm <- NA
                            }
                          } else stop("Precision issue with normalising constant")
                        } else {
                          ppas[[crit]]$PPAs <- exp(ppas[[crit]]$lPPA - log(norms))
                          ppas[[crit]]$norm <- norms
                        }
                        ppas["lPPA"] <- NULL
                        ppas
                      }, ppas = ppas, norms = norms, estimate = estimate, crit = crit)
                      names(ppas) <- names.priors
                      ppas
                    }, ppas = PPA_mat, norms = norms, names.priors = names(PPA_mat[[1]]), 
                      estimate = object$estimate, crit = 1)
                    PPAs.adjust <- PPAs.adjust[[1]]
                    # output results
                    crit <- 1
                    n_samp_mod <- ncol(PPAs.adjust[[1]][[crit]]$models_num)
                    models_num <- lapply(PPAs.adjust, function(x, crit) x[[crit]]$models_num, 
                      crit = crit)
                    hyps <- lapply(PPAs.adjust, function(x, crit) x[[crit]]$hyps, 
                      crit = crit)
                    ppas <- lapply(PPAs.adjust, function(x, crit) x[[crit]]$PPAs, 
                      crit = crit)
                    # sort into the correct order and amalgamate for printing
                    models_num <- lapply(models_num, function(models) {
                      # generate character representation of models
                      models <- apply(models, 1, function(mod, nsamp) numtochar(mod, 
                        nsamp), nsamp = ncol(models)/2)
                      models
                    })
                    if (length(models_num) > 1) {
                      # match up models
                      models <- do.call("c", models_num)
                      hyps <- do.call("c", hyps)
                      hyps <- hyps[!duplicated(models)]
                      models <- models[!duplicated(models)]
                      models_num <- lapply(models_num, function(mods, compmods) match(compmods, 
                        mods), compmods = models)
                      # now expand PPAs for each prior separately and fill in the blanks
                      ppas <- lapply(as.list(1:length(ppas)), function(i, ppas, orders) {
                        ppas <- ppas[[i]]
                        orders <- orders[[i]]
                        ppas <- ppas[orders[!is.na(orders)]]
                        orders <- match(orders, orders[!is.na(orders)])
                        ppas1 <- rep(NA, length(orders))
                        ppas1[orders[!is.na(orders)]] <- ppas
                        ppas1
                      }, ppas = ppas, orders = models_num)
                      models_num <- models
                    }
                    # append models to output data frame
                    output <- data.frame(models_num, ppas, hyps)
                    colnames(output) <- c("models", paste("(", object$priorPA, ")", 
                      sep = ""), "Null/Alt")
                  } else stop(paste("Criteria '", hyp_name[1], "' was not evaluated in this case", 
                    sep = ""))
                } else {
                  # check whether criteria is missing or not
                  if (!is.na(object$PPAs$nprior[3])) {
                    # now search for 'top' models
                    seqdata <- matrix(object$basedist[, ind], ncol = 1)
                    
                    # calculate models
                    nsamp <- nrow(seqdata)/4
                    if (nsamp == 1) 
                      stop("Need more than one sample")
                    
                    priorPA <- object$priorPA
                    
                    # produce prior specifications by generating (but not recording) each potential
                    # model sequentially
                    structure <- .Call("genmodels_priors", nsamp, PACKAGE = "seqmutprobs")
                    nstructure <- length(structure)
                    totmods <- structure[nstructure - 2]
                    nalt_ls <- structure[nstructure - 1]
                    nalt_s <- structure[nstructure]
                    nnull_ls <- totmods - nalt_ls
                    nnull_s <- totmods - nalt_s
                    structure <- structure[1:(nstructure - 3)]
                    # generate priorPAs for different hypotheses
                    lpriornull_ls <- (1 - priorPA)/nnull_ls
                    lprioralt_ls <- priorPA/nalt_ls
                    lpriornull_s <- (1 - priorPA)/nnull_s
                    lprioralt_s <- priorPA/nalt_s
                    lpriors <- cbind(lpriornull_ls, lprioralt_ls, lpriornull_s, lprioralt_s)
                    lpriors <- log(lpriors)
                    
                    if (hyp_name[1] == "less") 
                      nalt_s <- NA
                    if (hyp_name[1] == "stringent") 
                      nalt_ls <- NA
                    nprior <- c(nmod = totmods, nalt_less = nalt_ls, nalt_string = nalt_s)
                    # calculate log-threshold
                    logc <- log(c)
                    
                    # generate total number of columns required to store intermediate calcs
                    ntotcol <- sum(choose(nsamp, 1:nsamp))
                    
                    PPA_mat <- apply(seqdata, 2, function(seqs, nsamp, pstar, ntotcol, 
                      logc, lpriors, structure, criteria) {
                      # generate (10 x ncol)-matrix of intermediate values for calculating PPAs
                      lPDM_int_mat <- .C("calc_lPDM_int_fn", as.integer(nsamp), as.integer(seqs), 
                        as.double(pstar), lPDM_int_mat = as.double(numeric(10 * ntotcol)), 
                        PACKAGE = "seqmutprobs")$lPDM_int_mat
                      lPPAs <- apply(lpriors, 1, function(x, seqs, nsamp, pstar, 
                        ntotcol, logc, structure, criteria, lPDM_int_mat) {
                        # sort out and remove duplicates in 'lPDM_int_mat'
                        lPDM_int_mat1 <- lPDM_int_mat
                        lPDM_int_mat <- matrix(lPDM_int_mat, 10)
                        uni <- matrix(0, 10, nsamp * 10)
                        uni_ind <- matrix(0, 10, nsamp)
                        for (i in 1:nsamp) {
                          z <- unique(lPDM_int_mat[duplicated(lPDM_int_mat[, i]), 
                            i])
                          if (length(z) > 0) {
                            for (j in z) {
                              y <- which(lPDM_int_mat[, i] == j)
                              y1 <- y[1]
                              y <- y[2:length(y)]
                              lPDM_int_mat[y, i] <- (-1e+10)
                              uni[y1, (i - 1) * 10 + (1:length(y))] <- y - 1
                              uni_ind[y1, i] <- length(y)
                            }
                          }
                        }
                        lPDM_int_mat <- as.numeric(lPDM_int_mat)
                        uni <- as.numeric(uni)
                        uni_ind <- as.numeric(uni_ind)
                        if (criteria == "both" | criteria == "less") {
                          # now calculate PPAs according to approximation routine for LESS-STRINGENT
                          # criteria
                          lPPA_mat_ls <- .Call("calc_PPAs_approx_fn", nsamp, ntotcol, 
                            logc, x[1], x[2], lPDM_int_mat1, 0, structure, length(structure)/nsamp, 
                            1, uni, uni_ind, PACKAGE = "seqmutprobs")
                          totmods <- lPPA_mat_ls[length(lPPA_mat_ls)]
                          lPPA_mat_ls <- lPPA_mat_ls[1:(length(lPPA_mat_ls) - 1)]
                          models_num <- lPPA_mat_ls[1:(2 * nsamp * totmods)]
                          models_num <- matrix(models_num, nrow = totmods, byrow = T)
                          hyp <- lPPA_mat_ls[(2 * nsamp * totmods + 1):length(lPPA_mat_ls)]
                          lPPA_mat_ls <- hyp[(totmods + 1):length(hyp)]
                          multfact <- lPPA_mat_ls[(totmods + 1):length(lPPA_mat_ls)]
                          hyp <- hyp[1:totmods]
                          lPPA_mat_ls <- lPPA_mat_ls[1:totmods]
                          lPPA_mat_ls <- list(models_num = models_num, hyps = hyp, 
                            lPPA = lPPA_mat_ls)
                        }
                        if (criteria == "both" | criteria == "stringent") {
                          # now calculate PPAs according to approximation routine for LESS-STRINGENT
                          # criteria
                          lPPA_mat_s <- .Call("calc_PPAs_approx_fn", nsamp, ntotcol, 
                            logc, x[3], x[4], lPDM_int_mat1, 1, structure, length(structure)/nsamp, 
                            1, uni, uni_ind, PACKAGE = "seqmutprobs")
                          totmods <- lPPA_mat_s[length(lPPA_mat_s)]
                          lPPA_mat_s <- lPPA_mat_s[1:(length(lPPA_mat_s) - 1)]
                          models_num <- lPPA_mat_s[1:(2 * nsamp * totmods)]
                          models_num <- matrix(models_num, nrow = totmods, byrow = T)
                          hyp <- lPPA_mat_s[(2 * nsamp * totmods + 1):length(lPPA_mat_s)]
                          lPPA_mat_s <- hyp[(totmods + 1):length(hyp)]
                          multfact <- lPPA_mat_s[(totmods + 1):length(lPPA_mat_s)]
                          hyp <- hyp[1:totmods]
                          lPPA_mat_s <- lPPA_mat_s[1:totmods]
                          lPPA_mat_s <- list(models_num = models_num, hyps = hyp, 
                            lPPA = lPPA_mat_s)
                        }
                        if (criteria == "less") 
                          lPPA_mat_s <- NA
                        if (criteria == "stringent") 
                          lPPA_mat_ls <- NA
                        # output lists
                        list(less = lPPA_mat_ls, stringent = lPPA_mat_s)
                      }, seqs = seqs, nsamp = nsamp, pstar = pstar, ntotcol = ntotcol, 
                        logc = logc, structure = structure, criteria = criteria, 
                        lPDM_int_mat = lPDM_int_mat)
                      names(lPPAs) <- paste("priorPA_", priorPA, sep = "")
                      lPPAs
                    }, nsamp = nsamp, pstar = object$pstar, ntotcol = ntotcol, logc = logc, 
                      lpriors = lpriors, structure = structure, criteria = hyp_name[1])
                    # extract correct values from 'norms' if required
                    norms <- matrix(object$PPAs$norms[, ind], ncol = 1)
                    # calculate subset of top models with exact PPAs
                    PPAs.adjust <- lapply(as.list(1:ncol(norms)), function(i, ppas, 
                      norms, names.priors, estimate, crit) {
                      # extract correct sequence
                      ppas <- ppas[[i]]
                      norms <- matrix(norms[, i], nrow = 2)
                      # cycle through priors
                      ppas <- lapply(as.list(1:ncol(norms)), function(j, ppas, norms, 
                        estimate, crit) {
                        ppas <- ppas[[j]]
                        norms <- norms[crit, j]
                        # cycle through criteria
                        if (is.na(norms)) {
                          if (estimate == "top") {
                            ppa <- exp(ppas[[crit]]$lPPA)
                            if (length(which(!is.finite(log(ppa)))) > 0) {
                              norm <- ppa
                              prec <- 60
                              while (length(which(!is.finite(log(norm)))) > 0 & prec <= 
                                240) {
                                prec <- prec * 2
                                norm <- exp(mpfr(ppas[[crit]]$lPPA, prec))
                                if (length(which(!is.finite(log(norm)))) == 0) {
                                  norm <- sum(norm)
                                  ppa <- mpfr(ppas[[crit]]$lPPA, prec)
                                  ppa <- ppa - log(norm)
                                  ppa <- as.numeric(exp(ppa))
                                }
                              }
                              if (length(which(!is.finite(log(norm)))) > 0) 
                                stop("Precision issue with normalising constant")
                              ppas[[crit]]$PPAs <- ppa
                              ppas[[crit]]$norm <- NA
                            }
                          } else stop("Precision issue with normalising constant")
                        } else {
                          ppas[[crit]]$PPAs <- exp(ppas[[crit]]$lPPA - log(norms))
                          ppas[[crit]]$norm <- norms
                        }
                        ppas["lPPA"] <- NULL
                        ppas
                      }, ppas = ppas, norms = norms, estimate = estimate, crit = crit)
                      names(ppas) <- names.priors
                      ppas
                    }, ppas = PPA_mat, norms = norms, names.priors = names(PPA_mat[[1]]), 
                      estimate = object$estimate, crit = 2)
                    PPAs.adjust <- PPAs.adjust[[1]]
                    # output results
                    crit <- 2
                    n_samp_mod <- ncol(PPAs.adjust[[1]][[crit]]$models_num)
                    models_num <- lapply(PPAs.adjust, function(x, crit) x[[crit]]$models_num, 
                      crit = crit)
                    hyps <- lapply(PPAs.adjust, function(x, crit) x[[crit]]$hyps, 
                      crit = crit)
                    ppas <- lapply(PPAs.adjust, function(x, crit) x[[crit]]$PPAs, 
                      crit = crit)
                    # sort into the correct order and amalgamate for printing
                    models_num <- lapply(models_num, function(models) {
                      # generate character representation of models
                      models <- apply(models, 1, function(mod, nsamp) numtochar(mod, 
                        nsamp), nsamp = ncol(models)/2)
                      models
                    })
                    if (length(models_num) > 1) {
                      # match up models
                      models <- do.call("c", models_num)
                      hyps <- do.call("c", hyps)
                      hyps <- hyps[!duplicated(models)]
                      models <- models[!duplicated(models)]
                      models_num <- lapply(models_num, function(mods, compmods) match(compmods, 
                        mods), compmods = models)
                      # now expand PPAs for each prior separately and fill in the blanks
                      ppas <- lapply(as.list(1:length(ppas)), function(i, ppas, orders) {
                        ppas <- ppas[[i]]
                        orders <- orders[[i]]
                        ppas <- ppas[orders[!is.na(orders)]]
                        orders <- match(orders, orders[!is.na(orders)])
                        ppas1 <- rep(NA, length(orders))
                        ppas1[orders[!is.na(orders)]] <- ppas
                        ppas1
                      }, ppas = ppas, orders = models_num)
                      models_num <- models
                    }
                    # append models to output data frame
                    output <- data.frame(models_num, ppas, hyps)
                    colnames(output) <- c("models", paste("(", object$priorPA, ")", 
                      sep = ""), "Null/Alt")
                  } else stop(paste("Criteria '", hyp_name[1], "' was not evaluated in this case", 
                    sep = ""))
                }
            }
            # check whether the correct number of models can be output (in case using the
            # 'top' routine)
            if (num_mod > nrow(output)) 
                num_mod <- nrow(output)
            # sort into decreasing order according to PPA and extract top models
            output <- output[order(-output[, 2], output[, 1]), ]
            output <- output[1:num_mod, ]
            # print base distribution
            cat(paste("\nDistribution of bases at site ", i, ":\n", sep = ""))
            basedist <- matrix(object$basedist[, ind], 4)
            cons <- object$cons[i]
            bases <- c("a", "c", "g", "t")
            bases <- c(bases[bases != cons], bases[bases == cons])
            rownames(basedist) <- toupper(bases)
            colnames(basedist) <- object$samp_names
            basedistprop <- round(apply(basedist, 2, function(x) x/sum(x)), digits = 2)
            basedistprop.draw <- basedistprop
            basedistprop <- rbind(colnames(basedist), basedistprop)
            basedist <- rbind(colnames(basedist), basedist)
            write.table(format(basedist, justify = "centre"), quote = F, na = "", 
                col.names = F)
            # print base distribution as proportions (easier to visualise for NGS data)
            cat(paste("\nDistribution of bases at site", i, "as proportions:\n"))
            write.table(format(basedistprop, justify = "centre"), quote = F, na = "", 
                col.names = F)
            # print top num_mod models
            if (hyp_name[1] == "less") 
                hyp_name <- "LESS STRINGENT" else hyp_name <- "STRINGENT"
            cat(paste("\nTop ", num_mod, " models at site ", i, " based on '", hyp_name, 
                "' criterion:\n", sep = ""))
            output[, 2:ncol(output)] <- round(output[, 2:ncol(output)], digits = digits)
            # output1<-output
            output <- apply(output, 2, as.character)
            output <- rbind(colnames(output), output)
            output <- format(output, justify = "centre")
            write.table(output, quote = F, na = "", row.names = F, col.names = F)
            cat(paste("\nOverall PPA_SI for ", hyp_name, " criterion:\n", sep = ""))
            output <- sapply(object$hyp_PPAs, function(ppas, ind, crit) {
                if (crit == "LESS STRINGENT") 
                  x <- ppas[[1]][[ind]] else x <- ppas[[2]][[ind]]
                x
            }, ind = ind, crit = hyp_name)
            output <- round(output, digits = digits)
            output <- as.character(output)
            output <- rbind(sapply(strsplit(names(object$hyp_PPAs), "_"), function(x) paste(x[1], 
                " (", x[2], ")", sep = "")), output)
            output <- format(output, justify = "centre")
            write.table(output, quote = F, na = "", row.names = F, col.names = F)
            # add entropy calculations calculate normalised entropy measure
            dists <- object$basedist[, ind]
            dists <- matrix(dists, nrow = 4)
            entropy <- numeric(ncol(dists) - 1)
            # calculate entropy
            for (j in 2:ncol(dists)) {
                # produce normalised entropy
                s2 <- sum(dists[, j])
                s2 <- s2 * diag(4)
                ans1 <- apply(s2, 2, function(x, dists) entropy.fn(dists[, 1], x), 
                  dists = dists)
                entropy[j - 1] <- entropy.fn(dists[, 1], dists[, j])/max(ans1)
            }
            entropy <- c(entropy, mean(entropy), max(entropy))
            entropy <- round(entropy, digits = digits)
            entropy <- as.character(entropy)
            entropy <- matrix(entropy, nrow = 1)
            entropy <- rbind(c(paste(2:(ncol(entropy) - 1), ":1", sep = ""), "Mean", 
                "Max"), entropy)
            cat("\nK-L entropy measures\n")
            entropy <- format(entropy, justify = "centre")
            write.table(entropy, quote = F, na = "", row.names = F, col.names = F)
            
            # calculate shannon entropy
            shannon <- apply(dists, 2, function(x) {
                dists <- matrix(x, nrow = 4)
                # calculate normalised Shannon entropy
                ans <- apply(dists, 2, shannon.fn)
                ans <- ans/(-4 * 0.25 * log(0.25))
                ans
            })
            shannon <- c(shannon, mean(shannon), max(shannon))
            shannon <- round(shannon, digits = digits)
            shannon <- as.character(shannon)
            shannon <- matrix(shannon, nrow = 1)
            shannon <- rbind(c(paste(2:(ncol(shannon) - 1), ":1", sep = ""), "Mean", 
                "Max"), shannon)
            cat("\nShannon entropy measures\n")
            shannon <- format(shannon, justify = "centre")
            write.table(shannon, quote = F, na = "", row.names = F, col.names = F)
            
            # calculate K-L entropy relative to the prior
            klprior <- apply(dists, 2, function(dists) {
                # produce normalised entropy
                s <- sum(dists)
                s <- s * diag(4)
                ans1 <- apply(s, 2, klprior.fn)
                ans <- klprior.fn(dists)/max(ans1)
            })
            klprior <- 1 - klprior
            klprior <- c(klprior, mean(klprior), max(klprior))
            klprior <- round(klprior, digits = digits)
            klprior <- as.character(klprior)
            klprior <- matrix(klprior, nrow = 1)
            klprior <- rbind(c(paste(2:(ncol(klprior) - 1), ":1", sep = ""), "Mean", 
                "Max"), klprior)
            cat("\nK-L entropy measures (relative to prior)\n")
            klprior <- format(klprior, justify = "centre")
            write.table(klprior, quote = F, na = "", row.names = F, col.names = F)
            
            # draw if required
            if (draw == TRUE) {
                basedistprop.draw <- data.frame(Base = factor(rep(rownames(basedistprop.draw), 
                  times = ncol(basedistprop.draw)), levels = c("A", "C", "G", "T")), 
                  Sample = factor(rep(colnames(basedistprop.draw), each = nrow(basedistprop.draw)), 
                    levels = colnames(basedistprop.draw)), Proportion = matrix(basedistprop.draw, 
                    ncol = 1))
                temp.plot <- ggplot(data = basedistprop.draw, aes(Sample, Proportion, 
                  fill = Base)) + geom_bar(stat = "identity") + ggtitle(paste("Nucleotide site", 
                  site))
                print(temp.plot)
            }
        }
    }
    # if(!exists('output1')) output1<-NA return(output1)
} 
