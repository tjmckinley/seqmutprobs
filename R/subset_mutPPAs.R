#' Subset method for 'mutPPAs' objects
#' 
#' Subset method for 'mutPPAs' objects
#' 
#' @param x mutPPAs object
#' @param i a scalar indicating the nucleotide site to summarise
#' @param ... ignored
#' @rdname subset.mutPPAs
#' @method [ mutPPAs
#' @export [.mutPPAs

# subset method for 'mutPPAs' objects
"[.mutPPAs" <- function(x, i) {
    if (missing(i)) 
        return(x)
    if (is.null(i)) 
        return(x)
    subset <- i
    if (all(abs(subset) > 0 & abs(subset) < length(x$nucind))) {
        if (max(table(abs(subset))) > 1) 
            stop("Some identical negative and positive elements for subsetting")
        neg <- abs(subset[subset < 0])
        pos <- subset[subset > 0]
        # remove negative subsets
        if (length(neg) > 0) {
            x$nucind[neg] <- NA
            x$rem_sites <- c(x$rem_sites, neg)
            x$rem_sites <- sort(unique(x$rem_sites))
        }
        if (length(pos) > 0) {
            # retain positive sites
            temp <- which(!is.na(match(pos, x$rem_sites)))
            if (length(temp) > 0) {
                if (length(temp) == 1) 
                  temp.err <- paste("Site:", pos[temp], "already removed from analysis\n") else temp.err <- paste("Sites: ", paste(pos[temp[-length(temp)]], 
                  collapse = ", "), pos[temp[length(temp)]], " already removed from analysis\n", 
                  sep = "")
                pos <- pos[-temp]
            }
            if (length(pos) > 0) {
                x$nucind[-pos] <- NA
                x$rem_sites <- (1:length(x$nucind))[-pos]
                if (exists("temp.err")) 
                  cat(temp.err)
            } else stop(temp.err)
        }
    } else stop("Some elements of 'subset' were not in correct range")
    return(x)
}
 
