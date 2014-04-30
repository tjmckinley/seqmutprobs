# function to turn numerical representation of a single model to character form
# for printing
numtochar <- function(model, nsamp) {
    # extract models and model indicators
    inds <- model[(nsamp + 1):(2 * nsamp)]
    # convert any multiple M0s to correct indicator for printing
    mods <- model[1:nsamp]
    # inds[mods==0]<-inds[mods==0][1] convert to character
    mods_char <- as.character(mods)
    # generate counts from indictators
    counts <- table(inds)
    # produce character representation of models
    if (length(counts) != nsamp) {
        counts <- as.numeric(names(counts[counts > 1]))
        suffix <- letters[1:length(counts)]
        if (length(suffix[is.na(suffix)]) > 0) 
            stop("Too many combinations to produce character suffixes for models.")
        for (i in 1:length(counts)) mods_char[inds == counts[i]] <- paste(mods_char[inds == 
            counts[i]], suffix[i], sep = "")
        mods_char <- paste(mods_char, collapse = " ")
    } else mods_char <- paste(mods_char, collapse = " ")
    # return character representation
    mods_char
}
 
