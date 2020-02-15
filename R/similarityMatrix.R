#' Similarity score matrix
#' 
#' Calculates the similarity score of all pairwise comparison from a list of 
#' data frames.
#' 
#' @param productive.seqs A list data frames of productive sequences generated 
#' by the LymphoSeq function productiveSeq.  "count" and "aminoAcid" 
#' are a required columns.
#' @return A data frame of similarity scores calculated from all pairwise 
#' comparisons.  The similarity scores is a measure of the amount of overlap 
#' between two samples.  The value ranges from 0 to 1 where 1 indicates 
#' the sequence frequencies are identical in the two samples and 0 indicates 
#' no shared frequencies.
#' @seealso \code{\link{pairwisePlot}} for plotting results as a heat map.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' similarityMatrix(productive.seqs = productive.aa)
#' @export
similarityMatrix <- function(productive.seqs) {
    l <- length(productive.seqs)
    m <- matrix(nrow = l, ncol = l)
    rownames(m) <- names(productive.seqs)
    colnames(m) <- names(productive.seqs)
    for (i in 1:l) {
        for (j in 1:i) {
            m[i, j] <- similarityScore(productive.seqs[[i]], productive.seqs[[j]])
            m[j, i] <- m[i, j]
        }
    }
    m <- m[, colnames(m)[order(nchar(colnames(m)), colnames(m), decreasing = TRUE)]]
    m <- m[colnames(m)[order(nchar(colnames(m)), colnames(m), decreasing = TRUE)], ]
    d <- as.data.frame(m)
    return(d)
}
