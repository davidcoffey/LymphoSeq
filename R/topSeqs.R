#' Top sequences
#' 
#' Creates a data frame of a selected number of top productive sequences 
#' from a list of data frames.
#' 
#' @param productive.seqs A list data frames of productive sequences generated 
#' by the LymphoSeq function productiveSeq.  "frequencyCount" and "aminoAcid" 
#' are a required columns.
#' @param top The number of top productive sequences in each data frame to subset 
#' by their frequencies.
#' @return Returns a data frame of a selected number of top productive sequences 
#' from a list of data frames.
#' @seealso \code{\link{chordDiagramVDJ}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' top.seqs <- topSeqs(productive.seqs = productive.aa, top = 1)
#' @export
#' @importFrom plyr llply ldply
topSeqs <- function(productive.seqs, top = 1) {
    top.seqs <- plyr::llply(productive.seqs, function(x) 
        x[order(x$frequencyCount, decreasing = TRUE), ])
    top.seqs <- plyr::llply(productive.seqs, function(x) x[1:top, ])
    top.seqs <- plyr::ldply(top.seqs, data.frame)
    colnames(top.seqs)[1] <- "samples"
    return(top.seqs)
} 