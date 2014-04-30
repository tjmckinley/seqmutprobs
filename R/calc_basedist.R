# calculate distribution of bases and sort so that consensus base is in position
# 4
calc_basedist <- function(seqs, cons) {
    basedist <- apply(seqs, 2, function(x) {
        x <- x[!is.na(x)]
        c(length(x[x == 1]), length(x[x == 2]), length(x[x == 3]), length(x[x == 
            4]))
    })
    basedist <- rbind(basedist, cons)
    rows <- 1:4
    basedist <- apply(basedist, 2, function(x, rows) {
        cons <- x[length(x)]
        x <- x[1:(length(x) - 1)]
        x <- x[c(rows[-cons], cons)]
        x
    }, rows = rows)
    names(basedist) <- 1:ncol(basedist)
    basedist
}
 
