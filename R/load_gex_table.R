# Loads a SCAN-B gene expression matrix (samples in rows, genes in columns).
#
# When first called, the gene expression matrix in TSV format is read (slow!) and then saved
# as an R object in the same directory as the specified TSV file.  On subsequent calls the
# function tries to load the R object first, and when failing to do so falls back on reading
# the CSV file.
#
# @param gex.path path to a gene expression TSV file
# @return a matrix with genes in rows and samples in columns
# @export
load_gex_table <- function(gex.path) {
    path.components <- strsplit(gex.path, "\\.(?!.*\\.)", perl = TRUE)
    robj.path <- paste0(c(unlist(path.components)[1], "Rdata"), collapse=".")
    stopifnot(gex.path != robj.path)  # ensure we don't overwrite the original file
    if (file.exists(robj.path)) {
        writeLines(paste("Reading GEX cached Rdata file:", gex.path))
        gex <- get(load(robj.path))
    } else {
        writeLines(paste("Reading GEX text table:", gex.path))
        gex <- read.delim(gex.path, stringsAsFactors=F)
        rownames(gex) <- gex[, 1]
        gex <- gex[, -1]
        colnames(gex) <- gsub("^X(.+)", "\\1", colnames(gex))  # strip the leading X
        gex <- as.matrix(gex)
        writeLines(paste("Saving GEX table as Rdata cache file:", robj.path))
        save(gex, file=robj.path)
    }
    return(gex)
}