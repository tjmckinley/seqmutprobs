# function to convert alignment objects into pile-up tables
convert_align_pile <- function(seqdata, pstar = NULL, sites = NA, samp_names, genes, 
    reference) {
    #'seqdata' must be a list of 'alignment' objects from the 'seqinr' package
    if (!is.list(seqdata)) 
        stop("'seqdata' not a list of alignment objects")
    if (sum(sapply(lapply(seqdata, class), function(x) ifelse(x == "alignment", 0, 
        1))) > 0) 
        stop("'seqdata' not a list of alignment objects")
    
    nsamp <- length(seqdata)
    if (nsamp == 1) 
        stop("Need more than one sample")
    
    # extract number of sequences for summary
    nseq <- sapply(seqdata, function(x) x$nb)
    
    # now convert all alignments into character matrix format
    seqdata <- lapply(seqdata, function(x) {
        x <- as.matrix(x)
        x <- tolower(x)
    })
    # extract number of nucleotide sites for summary
    nnuc <- sapply(seqdata, function(x) ncol(x))
    
    # check that the number of nucleotide sites is the same across all samples
    if (length(unique(nnuc)) != 1) 
        stop("Different samples do not contain the same number of nucleotide sites")
    nnuc <- nnuc[1]
    
    # calculate consensus from initial sample
    cons <- calc_cons(seqdata[[1]], reference)
    # check for dodgy sites
    cons_tab <- names(table(cons))
    if (length(cons_tab[cons_tab != "a" & cons_tab != "c" & cons_tab != "g" & cons_tab != 
        "t" & cons_tab != "-"]) > 0) 
        stop("A nucleotide site exists in the initial sample that does not have a valid consensus (e.g. none of 'a', 'c', 'g', 't' or '-')")
    
    # check for and remove insertions if necessary
    ins_loc <- grep("-", cons)
    full_cons <- cons
    if (length(ins_loc) > 0) {
        # cat(paste('Positions:',paste(ins_loc,collapse=' '),'removed as possible
        # insertion sites (i.e. initial sample consensus of '-')\n'))
        seqdata <- lapply(seqdata, function(x, ins_loc) x[, -ins_loc], ins_loc = ins_loc)
        cons <- cons[-ins_loc]
    } else ins_loc <- NA
    
    # now convert to numeric format (s.t. {acgt}->{1234}) for quicker processing
    seqdata <- lapply(seqdata, function(x) t(apply(x, 1, s2n)) + 1)
    cons <- s2n(cons) + 1
    
    # now calculate distribution of bases at each site and bind together
    seqdata <- do.call("rbind", lapply(seqdata, function(x, cons) calc_basedist(x, 
        cons), cons = cons))
    
    # if 'pstar' not specified then calculate
    if (is.null(pstar)) 
        pstar <- sum(seqdata[(1:nrow(seqdata))[(1:nrow(seqdata))%%4 != 0], ])/sum(seqdata)
    
    # now remove sites that aren't selected for testing
    if (!is.na(sites[1])) {
        rem_sites <- (1:length(full_cons))[-sites]
        rems <- (as.numeric(colnames(seqdata)) %in% rem_sites)
        rems <- which(rems == T)
        seqdata <- seqdata[, -rems]
        if (is.null(ncol(seqdata))) 
            seqdata <- matrix(seqdata, ncol = 1)
        cons <- cons[-rems]
        if (!is.na(ins_loc[1])) {
            rem_sites <- c(rem_sites, ins_loc)
            rem_sites <- unique(rem_sites)
        }
    } else {
        if (!is.na(ins_loc[1])) 
            rem_sites <- ins_loc else rem_sites <- NA
    }
    
    # remove duplicate columns and create indicator to reduce memory requirements
    seqdata_temp <- seqdata
    seqdata <- seqdata[, !duplicated(seqdata, MARGIN = 2)]
    if (is.null(ncol(seqdata))) 
        seqdata <- matrix(seqdata, ncol = 1)
    colnames(seqdata) <- 1:ncol(seqdata)
    seqind_temp <- apply(seqdata_temp, 2, function(x, y) (1:ncol(y))[apply(y, 2, 
        function(y, x) all(x == y), x = x)], y = seqdata)
    # add back in insertion information so that 'seqind' is the same length as the
    # original sequence
    seqind <- numeric(nnuc)
    if (!is.na(rem_sites[1])) {
        seqind[rem_sites] <- NA
        seqind[-rem_sites] <- seqind_temp
    } else seqind <- seqind_temp
    rm(seqdata_temp, seqind_temp)
    
    # output data
    seqs <- list(seqdata = seqdata, seqind = seqind, pstar = pstar, rem_sites = rem_sites, 
        full_cons = full_cons, nnuc = nnuc, nseq = nseq, nsamp = nsamp, samp_names = samp_names, 
        sites = sites, genes = genes)
}
 
