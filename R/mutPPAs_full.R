# function for calculating PPAs from output of 'mutlPDMs_full'
mutPPAs_full <- function(mutlPDM, priorPA = c(0.001, 0.01, 0.05), criteria = c("both", 
    "stringent", "less"), mc.cores = 1, ...) {
    # sort priorPA (important for summary)
    priorPA <- sort(priorPA)
    
    # generate null and alternative hypotheses based on model structure and criteria
    if (criteria[1] == "both") {
        hyp_lessstring <- .C("calc_lessstring_fn", as.integer(length(mutlPDM$nseq)), 
            as.integer(nrow(mutlPDM$models_num)), as.integer(mutlPDM$models_num), 
            hyp = as.integer(numeric(nrow(mutlPDM$models_num))), PACKAGE = "seqmutprobs")$hyp
        hyp_string <- .C("calc_string_fn", as.integer(length(mutlPDM$nseq)), as.integer(nrow(mutlPDM$models_num)), 
            as.integer(mutlPDM$models_num), hyp = as.integer(numeric(nrow(mutlPDM$models_num))), 
            PACKAGE = "seqmutprobs")$hyp
    } else {
        if (criteria[1] == "stringent") {
            hyp_lessstring <- NA
            hyp_string <- .C("calc_string_fn", as.integer(length(mutlPDM$nseq)), 
                as.integer(nrow(mutlPDM$models_num)), as.integer(mutlPDM$models_num), 
                hyp = as.integer(numeric(nrow(mutlPDM$models_num))), PACKAGE = "seqmutprobs")$hyp
        } else {
            hyp_lessstring <- .C("calc_lessstring_fn", as.integer(length(mutlPDM$nseq)), 
                as.integer(nrow(mutlPDM$models_num)), as.integer(mutlPDM$models_num), 
                hyp = as.integer(numeric(nrow(mutlPDM$models_num))), PACKAGE = "seqmutprobs")$hyp
            hyp_string <- NA
        }
    }
    # produce priorPAs for different criteria
    if (criteria[1] == "both" | criteria[1] == "less") {
        # calculate how priorPAs are split across models for less stringent hypothesis
        nalt <- sum(hyp_lessstring)
        nnull <- length(hyp_lessstring) - nalt
        p_lessstring <- apply(matrix(priorPA, 1), 2, function(p, nalt, nnull, hyp) {
            palt <- p/nalt
            pnull <- (1 - p)/nnull
            probs <- hyp
            probs[hyp == 1] <- palt
            probs[hyp == 0] <- pnull
            probs
        }, nalt = nalt, nnull = nnull, hyp = hyp_lessstring)
        p_lessstring <- log(p_lessstring)
        p_lessstring <- matrix(p_lessstring, ncol = length(priorPA))
        # calculate PPAs for less stringent hypothesis
        PPAs_lessstring <- mclapply(as.list(1:ncol(p_lessstring)), function(i, pPA, 
            lPDM) {
            pPA <- pPA[, i]
            PPA <- apply(lPDM, 2, function(l, p) {
                l <- l + p
                ppa <- exp(l)
                if (length(which(!is.finite(log(ppa)))) > 0) {
                  prec <- 60
                  while (length(which(!is.finite(log(ppa)))) > 0 & prec <= 240) {
                    prec <- prec * 2
                    ppa <- mpfr(l, prec)
                    ppa <- exp(ppa)
                    if (length(which(!is.finite(log(ppa)))) == 0) {
                      norm <- sum(ppa)
                      ppa <- as.numeric(exp(log(ppa) - log(norm)))
                      return(list(ppa = ppa, norm = NA))
                    }
                  }
                  stop("Precision issue with normalising constant")
                } else {
                  norm <- sum(ppa)
                  ppa <- exp(l - log(norm))
                  return(list(ppa = ppa, norm = norm))
                }
            }, p = pPA)
            norm <- do.call("c", lapply(PPA, function(ppa) ppa$norm))
            PPA <- do.call("cbind", lapply(PPA, function(ppa) ppa$ppa))
            PPA <- list(norm = norm, PPAs = PPA)
            PPA
        }, pPA = p_lessstring, lPDM = mutlPDM$lPDM, mc.cores = mc.cores)
        # retain intermediate prior calculations to include in 'mutPPAs' class
        nalt_less <- nalt
    }
    if (criteria[1] == "both" | criteria[1] == "stringent") {
        # calculate how priorPAs are split across models for stringent hypothesis
        nalt <- sum(hyp_string)
        nnull <- length(hyp_string) - nalt
        p_string <- apply(matrix(priorPA, 1), 2, function(p, nalt, nnull, hyp) {
            palt <- p/nalt
            pnull <- (1 - p)/nnull
            probs <- hyp
            probs[hyp == 1] <- palt
            probs[hyp == 0] <- pnull
            probs
        }, nalt = nalt, nnull = nnull, hyp = hyp_string)
        p_string <- log(p_string)
        p_string <- matrix(p_string, ncol = length(priorPA))
        # calculate PPAs for stringent hypothesis
        PPAs_string <- mclapply(as.list(1:ncol(p_string)), function(i, pPA, lPDM) {
            pPA <- pPA[, i]
            PPA <- apply(lPDM, 2, function(l, p) {
                l <- l + p
                ppa <- exp(l)
                if (length(which(!is.finite(log(ppa)))) > 0) {
                  prec <- 60
                  while (length(which(!is.finite(log(ppa)))) > 0 & prec <= 240) {
                    prec <- prec * 2
                    ppa <- mpfr(l, prec)
                    ppa <- exp(ppa)
                    if (length(which(!is.finite(log(ppa)))) == 0) {
                      norm <- sum(ppa)
                      ppa <- as.numeric(exp(log(ppa) - log(norm)))
                      return(list(ppa = ppa, norm = NA))
                    }
                  }
                  stop("Precision issue with normalising constant")
                } else {
                  norm <- sum(ppa)
                  ppa <- exp(l - log(norm))
                  return(list(ppa = ppa, norm = norm))
                }
            }, p = pPA)
            norm <- do.call("c", lapply(PPA, function(ppa) ppa$norm))
            PPA <- do.call("cbind", lapply(PPA, function(ppa) ppa$ppa))
            PPA <- list(norm = norm, PPAs = PPA)
            PPA
        }, pPA = p_string, lPDM = mutlPDM$lPDM, mc.cores = mc.cores)
        # retain intermediate prior calculations to include in 'mutPPAs' class
        nalt_string <- nalt
    }
    if (criteria[1] == "less") 
        nalt_string <- NA
    if (criteria[1] == "stringent") 
        nalt_less <- NA
    # condense into correct format
    if (criteria[1] == "both") 
        PPAs <- lapply(as.list(1:length(PPAs_lessstring)), function(i, ppas_ls, ppas_s) list(less = ppas_ls[[i]], 
            stringent = ppas_s[[i]]), ppas_ls = PPAs_lessstring, ppas_s = PPAs_string) else {
        if (criteria[1] == "less") 
            PPAs <- lapply(as.list(1:length(PPAs_lessstring)), function(i, ppas_ls) list(less = ppas_ls[[i]], 
                stringent = list(norm = NA, PPAs = NA)), ppas_ls = PPAs_lessstring) else PPAs <- lapply(as.list(1:length(PPAs_string)), function(i, ppas_s) list(less = list(norm = NA, 
            PPAs = NA), stringent = ppas_s[[i]]), ppas_s = PPAs_string)
    }
    # output information
    norms <- lapply(PPAs, function(ppas) {
        norms <- lapply(ppas, function(ppas) ppas$norm)
        norms
    })
    names(norms) <- paste("priorPA_", priorPA, sep = "")
    PPAs <- lapply(PPAs, function(ppas) {
        PPAs <- lapply(ppas, function(ppas) ppas$PPAs)
        PPAs
    })
    names(PPAs) <- paste("priorPA_", priorPA, sep = "")
    # re-arrange to match output from 'mutPPAs_top'
    PPAs <- lapply(as.list(1:ncol(mutlPDM$basedist)), function(i, ppas) {
        lapply(ppas, function(ppas, i) {
            lapply(ppas, function(ppas, i) {
                if (!is.na(ppas[1])) 
                  ppas <- ppas[, i] else ppas <- NA
                ppas
            }, i = i)
        }, i = i)
    }, ppas = PPAs)
    names(PPAs) <- colnames(mutlPDM$basedist)
    # re-arrange to match output from 'mutPPAs_top'
    norms <- lapply(as.list(1:ncol(mutlPDM$basedist)), function(i, norms) {
        lapply(norms, function(norms, i) {
            lapply(norms, function(norms, i) {
                if (!is.na(norms[1])) 
                  norms <- as.numeric(norms[i]) else norms <- NA
                norms
            }, i = i)
        }, i = i)
    }, norms = norms)
    names(norms) <- colnames(mutlPDM$basedist)
    hyps <- list(hyp_lessstring, hyp_string)
    nmods <- sapply(hyps, length)
    if (criteria[1] == "less") 
        nmods[2] <- NA
    if (criteria[1] == "stringent") 
        nmods[1] <- NA
    nprior <- c(nmod = nmods[!is.na(nmods)][1], nalt_less = nalt_less, nalt_string = nalt_string)
    PPAs <- list(models_num = mutlPDM$models_num, hyps = hyps, nprior = nprior, PPAs = PPAs, 
        norms = norms)
    names(PPAs$hyps) <- c("less", "stringent")
    
    # reset names of sequences
    names(PPAs$PPAs) <- as.character(1:length(PPAs$PPAs))
    names(PPAs$norms) <- as.character(1:length(PPAs$norms))
    
    # output object of class 'mutPPAs'
    output <- mutlPDM
    output["lPDM"] <- NULL
    output["models_num"] <- NULL
    output$PPAs <- PPAs
    output$priorPA <- priorPA
    output$estimate <- "full"
    output$supp_output <- FALSE
    output <- output[c(1:8, 10, 9, 11:length(output))]
    
    # return 'mutPPAs' object
    output
} 
