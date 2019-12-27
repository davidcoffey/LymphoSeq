#' Bhattacharyya coefficient
#' 
#' Calculates the Bhattacharyya coefficient of two samples.
#' 
#' @param sample1 A data frame consisting of frequencies of antigen receptor 
#' sequences.  "frequencyCount" is a required column.
#' @param sample2 A data frame consisting of frequencies of antigen receptor 
#' sequences.  "frequencyCount" is a required column.
#' @return Returns the Bhattacharyya coefficient, a measure of the amount of 
#' overlap between two samples.  The value ranges from 0 to 1 where 1 indicates 
#' the sequence frequencies are identical in the two samples and 0 indicates no 
#' shared frequencies. 
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list, aggregate = "aminoAcid")
#' 
#' bhattacharyyaCoefficient(productive.aa[["TRB_Unsorted_32"]], 
#'    productive.aa[["TRB_Unsorted_83"]])
#' @seealso \code{\link{bhattacharyyaMatrix}}
#' @export
bhattacharyyaCoefficient <- function(sample1, sample2) {
    p <- sample1[which(sample1[, "aminoAcid"] %in% sample2[, "aminoAcid"]), ]
    p$frequencyCount <- p$frequencyCount/100
    q <- sample2[which(sample2[, "aminoAcid"] %in% sample1[, "aminoAcid"]), ]
    q$frequencyCount <- q$frequencyCount/100
    p <- p[match(q$aminoAcid, p$aminoAcid),]
    s <- p[, "frequencyCount"] * q[, "frequencyCount"]
    bc <- sum(sqrt(s))
    return(bc)
}