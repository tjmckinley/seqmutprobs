# function to create pileup tables from data files
create.pileup <- function(filenames, format = c("mase", "clustal", "phylip", "fasta", 
    "msf", "bam", "pileup"), pstar = NULL, genes = NA, sites = NA, samp_names = NULL, 
    cov_thresh = 5, reference = NULL, ref_file, mc.cores = 1, ...) {
    #'filenames' is a list of characters each containing paths to files we wish to import
    #'format' describes the formats of the files we want to import (all files must be the same format)
    #'pstar' is overall mutation probability 
    #'genes' specifies which segments to extract from any BAM files (only relevant for BAM files)
    #'sites' specify the sites to extract from the alignments
    #'samp_names' is vector of sample names
    #'cov_thresh' is the coverage threshold for extracting sites from NGS data
    #'reference' is a either a single character vector or a list containing the reference sequences
    # for different Gene segments if necessary
    #'ref_file' is path to reference file
    #'mc.cores' is the number of cores to use for parallel processing
    
    # read in sequences
    if (format[1] != "bam" & format[1] != "pileup") {
        if (!is.null(pstar[1]) & length(pstar) > 1) 
            stop("'pstar' specified incorrectly (must be of length 1 or NULL")
        if (!is.character(reference)) 
            stop("'reference' in wrong format")
        # read in a list of alignment objects
        seqs <- lapply(filenames, function(files, format) read.alignment(files, format = format), 
            format = format[1])
        # generate sample names
        if (is.null(samp_names[1])) {
            names(seqs) <- unlist(filenames)
            samp_names <- names(seqs)
        } else {
            if (length(samp_names) != length(seqs)) 
                stop("'samp_names' do not match 'filenames'") else names(seqs) <- samp_names
        }
        seqnames <- lapply(seqs, names)
        names(seqnames) <- samp_names
        # generate Gene segment names
        if (is.na(genes[1])) 
            genes <- "Segment 1" else if (length(genes) != 1) 
            stop("'genes' argument not of length 1")
        # convert list of 'alignment' objects to pile-up tables
        seqs <- convert_align_pile(seqs, pstar, sites, samp_names, genes, reference)
    } else {
        if (is.na(cov_thresh[1])) 
            stop("No coverage threshold given")
        # generate pile-up tables using 'samtools'
        if (format[1] == "bam") {
            if (!file.exists(paste(ref_file, ".fai", sep = ""))) 
                system(paste("samtools faidx", ref_file))
            pathtoperl <- paste(system.file(package = "seqmutprobs"), "/Perl/pileup2csv.pl", 
                sep = "")
            seqs <- mclapply(filenames, function(files, ref) read.table(pipe(paste("samtools mpileup -BQ0 -d10000000 -f ", 
                ref, " ", files, " | perl ", pathtoperl, sep = "")), header = TRUE, 
                quote = "\"", sep = "\t"), ref = ref_file, mc.cores = mc.cores)
        } else {
            pathtoperl <- paste(system.file(package = "seqmutprobs"), "/Perl/pileup2csv.pl", 
                sep = "")
            seqs <- lapply(filenames, function(files) read.table(pipe(paste("perl", 
                pathtoperl, files, sep = " ")), header = TRUE, quote = "\"", sep = "\t"))
        }
        
        if (is.null(samp_names[1])) {
            names(seqs) <- unlist(filenames)
            samp_names <- names(seqs)
        } else {
            if (length(samp_names) != length(seqs)) 
                stop("'samp_names' do not match 'filenames'") else names(seqs) <- samp_names
        }
        seqnames <- lapply(seqs, function(x) as.character(unique(x$NAME)))
        names(seqnames) <- samp_names
        
        # extract relevant genes from BAM files
        if (is.na(genes[1])) 
            genes <- seqnames[[1]]
        # check 'pstar' is either NULL, length 1 (in which case applied to all segments)
        # or of the same length as 'genes'
        if (!is.null(pstar[1]) & length(pstar) > 1 & length(pstar) != length(genes)) 
            stop("'pstar' is not the correct length (i.e. NULL, length 1 or the same length as 'genes'")
        
        matches <- lapply(seqnames, function(nam, genes) match(genes, nam), genes = genes)
        # remove any Gene segments not present in initial sample
        matchesini <- genes[!is.na(matches[[1]])]
        if (!is.null(pstar[1]) & length(pstar) > 1) 
            pstar <- pstar[!is.na(matches[[1]])]
        mismatches <- genes[is.na(matches[[1]])]
        # exit if no matches found in initial sample
        if (length(matchesini) < 1) 
            stop("No matching segments in initial sample")
        # otherwise carry on
        if (length(mismatches) > 0) {
            if (length(mismatches) > 2) 
                cat(paste(paste(paste(mismatches[1:(length(mismatches) - 1)], collapse = ", "), 
                  mismatches[length(mismatches)], collapse = " and "), "not present in initial sample and so are removed\n")) else cat(paste(paste(mismatches, collapse = "and"), "are only present in initial sample and so are removed\n"))
        }
        matches <- lapply(matches, function(x, matchesini) x[!is.na(matchesini)], 
            matchesini = matches[[1]])
        seqnames <- lapply(as.list(1:length(seqnames)), function(i, nam, matches) {
            nam <- nam[[i]]
            matches <- matches[[i]]
            nam[matches[!is.na(matches)]]
        }, nam = seqnames, matches = matches)
        # check at least two samples for each segment
        matches <- do.call("c", seqnames)
        matches <- table(matches)
        mismatches <- matches[matches < 2]
        if (!is.null(pstar[1]) & length(pstar) > 1) 
            pstar <- pstar[matches[matches < 2]]
        matches <- matches[matches >= 2]
        if (length(matches) < 1) 
            stop("No segments have more than one sample")
        if (length(mismatches) > 0) {
            if (length(mismatches) > 2) 
                cat(paste(paste(paste(names(mismatches[1:(length(mismatches) - 1)]), 
                  collapse = ", "), names(mismatches[length(mismatches)]), collapse = " and "), 
                  "are present in less than two samples and so are removed\n")) else cat(paste(paste(names(mismatches), collapse = "and"), "are present in less than two samples and so are removed\n"))
        }
        seqnames <- lapply(seqnames, function(nam, mismatch) {
            for (i in mismatch) nam <- nam[nam != i]
            nam
        }, mismatch = names(mismatches))
        names(seqnames) <- names(seqs)
        seqdata <- seqs
        rm(seqs)
        
        # remove sites with a low coverage threshold
        seqdata <- lapply(seqdata, function(x, cov_thresh) x[apply(x[, 5:8], 1, sum) > 
            cov_thresh, ], cov_thresh = cov_thresh)
        
        # merge pile-up tables across individuals for each gene segment
        genes <- seqnames[[1]]
        seqs <- list(length(genes))
        ref_match <- match(genes, names(reference))
        if (length(ref_match[is.na(ref_match)]) > 0) 
            stop("Reference sequence not available for some gene segments")
        # extract and sort reference sequences into correct order
        reference <- reference[ref_match]
        # now extract information by Gene segment
        for (i in 1:length(genes)) {
            # extract correct segment for each sample
            seqs[[i]] <- lapply(seqdata, function(seqs, nam) seqs[seqs$NAME == nam, 
                ], nam = genes[i])
            # remove any missing samples
            seqs[[i]] <- seqs[[i]][sapply(seqs[[i]], length) > 0]
            # figure out consensus and merge pileup tables together
            seqs[[i]] <- lapply(seqs[[i]], function(seqs, nref) {
                output <- matrix(NA, 4, nref)
                output[, seqs$POS] <- t(seqs[, 5:8])
                output
            }, nref = length(reference[[i]]))
            cons <- seqs[[i]][[1]]
            cons <- sapply(as.list(1:ncol(cons)), function(i, x) {
                ref <- which(c("a", "c", "g", "t") == x[5, i])
                x <- as.numeric(x[1:4, i])
                if (!is.na(sum(x))) {
                  y <- which(x == max(x))
                  if (length(y) > 1) {
                    if (length(!is.na(match(y, ref))) == 1) 
                      cons <- ref else {
                      cat(paste("Can't generate a consensus for site ", i, ",\nand can't match to reference in a sensible manner,\nso this site has been removed:\nREF = ", 
                        c("A", "C", "G", "T")[ref], "\n", sep = ""))
                      cons <- NA
                    }
                  } else cons <- y
                } else cons <- NA
                cons
            }, x = rbind(cons, reference[[i]]))
            
            # calculate length of nucleotide sequence and number of samples
            nnuc <- length(reference[[i]])
            nsamp <- length(seqs[[i]])
            
            # merge samples together
            seqs[[i]] <- do.call("rbind", seqs[[i]])
            seq_nam <- 1:ncol(seqs[[i]])
            
            # check for and remove invalid sites
            ins_loc <- which(is.na(apply(seqs[[i]], 2, sum)) | is.na(cons))
            full_cons <- cons
            full_cons <- n2s(full_cons - 1)
            if (length(ins_loc) > 0) {
                seqs[[i]] <- seqs[[i]][, -ins_loc]
                seq_nam <- seq_nam[-ins_loc]
                cons <- cons[-ins_loc]
            } else ins_loc <- NA
            
            if (length(cons) == 0) 
                seqs[[i]] <- NA else {
                # record range of number of sequences
                nseq <- list(nsamp)
                for (j in 1:nsamp) nseq[[j]] <- range(apply(seqs[[i]][1:4 + (j - 
                  1) * 4, ], 2, sum))
                
                # reorder rows so that consensus base is in row 4
                seqs[[i]] <- sapply(1:ncol(seqs[[i]]), function(i, x) {
                  x <- x[, i]
                  cons <- x[length(x)]
                  x <- x[-length(x)]
                  x <- matrix(x, 4)
                  x <- rbind(x[-cons, ], x[cons, ])
                  x <- as.numeric(x)
                  x
                }, x = rbind(seqs[[i]], cons))
                
                # if 'pstar' not specified then calculate
                if (is.null(pstar[1])) 
                  pstar1 <- sum(seqs[[i]][(1:nrow(seqs[[i]]))[(1:nrow(seqs[[i]]))%%4 != 
                    0], ])/sum(seqs[[i]]) else {
                  if (length(pstar) == 1) 
                    pstar1 <- pstar else pstar1 <- pstar[i]
                }
                # now remove sites that aren't selected for testing
                if (!is.na(sites[1])) {
                  rem_sites <- (1:length(full_cons))[-sites]
                  rems <- (seq_nam %in% rem_sites)
                  rems <- which(rems == T)
                  seq_nam <- seq_nam[-rems]
                  seqs[[i]] <- seqs[[i]][, -rems]
                  if (is.null(ncol(seqs[[i]]))) 
                    seqs[[i]] <- matrix(seqs[[i]], ncol = 1)
                  cons <- cons[-rems]
                  if (!is.na(ins_loc[1])) {
                    rem_sites <- c(rem_sites, ins_loc)
                    rem_sites <- unique(rem_sites)
                  }
                  colnames(seqs[[i]]) <- seq_nam
                } else {
                  if (!is.na(ins_loc[1])) 
                    rem_sites <- ins_loc else rem_sites <- NA
                }
                
                # remove duplicate columns and create indicator to reduce memory requirements
                seqdata_temp <- seqs[[i]]
                seqs[[i]] <- seqs[[i]][, !duplicated(seqs[[i]], MARGIN = 2)]
                if (is.null(ncol(seqs[[i]]))) 
                  seqs[[i]] <- matrix(seqs[[i]], ncol = 1)
                colnames(seqs[[i]]) <- 1:ncol(seqs[[i]])
                seqind_temp <- apply(seqdata_temp, 2, function(x, y) (1:ncol(y))[apply(y, 
                  2, function(y, x) all(x == y), x = x)], y = seqs[[i]])
                # add back in insertion information so that 'seqind' is the same length as the
                # original sequence
                seqind <- numeric(nnuc)
                if (!is.na(rem_sites[1])) {
                  seqind[rem_sites] <- NA
                  seqind[-rem_sites] <- seqind_temp
                } else seqind <- seqind_temp
                rm(seqdata_temp, seqind_temp)
                
                seqs[[i]] <- list(seqdata = seqs[[i]], seqind = seqind, pstar = pstar1, 
                  rem_sites = rem_sites, full_cons = full_cons, nnuc = nnuc, nseq = nseq, 
                  nsamp = nsamp, sites = sites, samp_names = samp_names, genes = genes[i])
            }
        }
        names(seqs) <- genes
    }
    seqs
} 
