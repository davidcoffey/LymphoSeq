#' Common sequences in two or more samples
#' 
#' Creates a data frame of the common sequences in two or more samples, reporting 
#' their frequencies in each.
#' 
#' @param samples A character vector of two or more sample names in 
#' productive.aa.
#' @param productive.aa A list of productive amino acid sequences generated
#' by the LymphoSeq function productiveSeq where aggregate = "aminoAcid".
#' @return Returns a data frame of the common sequences between two or more files 
#' displaying their frequencies in each.
#' @seealso \code{\link{commonSeqsVenn}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' commonSeqs(samples = c("TRB_Unsorted_0", "TRB_Unsorted_32"), 
#'    productive.aa = productive.aa)
#' @export
#' @importFrom plyr llply
#' @import ggplot2
commonSeqs <- function(samples, productive.aa) {
    aminoAcid <- plyr::llply(productive.aa[samples], 
                             function(x) x[, "aminoAcid"])
    common <- Reduce(intersect, aminoAcid)
    seq.matrix <- seqMatrix(productive.aa[samples], common)
    seq.matrix$numberSamples <- NULL
    return(seq.matrix)
}