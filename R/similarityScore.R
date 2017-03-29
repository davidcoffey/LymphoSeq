#' Similarity score
#' 
#' Calculates the similarity score of two samples.
#' 
#' @param sample1 A data frame consisting of frequencies of antigen receptor 
#' sequences.  "aminoAcid" and "count" are a required columns.
#' @param sample2 A data frame consisting of frequencies of antigen receptor 
#' sequences.  "aminoAcid" and "count" are a required columns.
#' @return Returns the similarity score, a measure of the amount of 
#' overlap between two samples.  The value ranges from 0 to 1 where 1 indicates 
#' the sequence frequencies are identical in the two samples and 0 
#' indicates no shared frequencies.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list, aggregate = "aminoAcid")
#' 
#' similarityScore(productive.aa[["TRB_Unsorted_32"]], productive.aa[["TRB_Unsorted_83"]])
#' @seealso \code{\link{similarityMatrix}}
#' @export
similarityScore <- function(sample1, sample2) {
    s1 <- sum(sample1[which(sample1[, "aminoAcid"] %in% sample2[, "aminoAcid"]), 
                      "count"])
    s2 <- sum(sample2[which(sample2[, "aminoAcid"] %in% sample1[, "aminoAcid"]), 
                      "count"])
    return((s1 + s2)/(sum(sample1[, "count"]) + sum(sample2[, "count"])))
}