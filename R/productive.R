#' Productive sequences
#' 
#' Remove unproductive CDR3 sequences from a single data frame.
#' 
#' @param sample A data frame consisting of antigen receptor sequencing data.  
#' "aminoAcid", "count", and "frequencyCount" are required columns.
#' @param aggregate Indicates whether the values of "count", "frequencyCount", 
#' and "esimatedNumberGenomes" should be aggregated by amino acid or nucleotide 
#' sequence.  Acceptable values are "aminoAcid" or "nucleotide".  If "aminoAcid" 
#' is selected, then the resulting data frame will have columns corresponding to 
#' "aminoAcid", "count", "frequnecyCount", and "estimatedNumberGenomes" (if this 
#' column is available).  If "nucleotide" is selected then all columns in the 
#' original data frame will be present in the outputted data frame.  The difference in 
#' output is due to the fact that the same amino acid CDR3 sequence may be 
#' encoded by multiple unique nucleotide sequences with differing V, D, and J 
#' genes. 
#' @return Returns a data frame of productive amino acid sequences with recomputed 
#' values for "count", "frequencyCount", and "esimatedNumberGenomes".  A 
#' productive sequences is defined as a sequence that is in frame and does not 
#' have an early stop codon.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive <- productive(sample = file.list[["TCRB_Day32_Unsorted"]], aggregate = "aminoAcid")
#' @seealso \code{\link{productiveSeq}}
#' @export
#' @importFrom stats aggregate
productive <- function(sample, aggregate = "aminoAcid") {
    productive.seqs <- sample[which(sample[, "aminoAcid"] != "" &
                                        !grepl("\\*", sample[, "aminoAcid"])), ]
    if(nrow(productive.seqs) != 0){
        if (any(grepl("estimatedNumberGenomes", colnames(productive.seqs)))) {
            if (aggregate == "aminoAcid") {
                count <- aggregate(count ~ aminoAcid, data = productive.seqs, FUN = sum)
                productive.seqs$estimatedNumberGenomes[is.na(as.integer(productive.seqs$estimatedNumberGenomes))] <- 0
                estimatedNumberGenomes <- aggregate(estimatedNumberGenomes ~ aminoAcid, 
                                                    data = productive.seqs, FUN = sum)
                merged <- merge(count, estimatedNumberGenomes, by = "aminoAcid")
                merged$frequencyCount <- merged$count/sum(count$count) * 100
            }
            if (aggregate == "nucleotide") {
                count <- aggregate(count ~ nucleotide, data = productive.seqs, FUN = sum)
                productive.seqs$estimatedNumberGenomes[is.na(as.integer(productive.seqs$estimatedNumberGenomes))] <- 0
                estimatedNumberGenomes <- aggregate(estimatedNumberGenomes ~ nucleotide, 
                                                    data = productive.seqs, FUN = sum)
                productive.seqs$count <- NULL
                productive.seqs$estimatedNumberGenomes <- NULL
                merged <- Reduce(function(x, y) merge(x, y, by = "nucleotide"), 
                                 list(productive.seqs, count, estimatedNumberGenomes))
                merged$frequencyCount <- merged$count/sum(count$count) * 100
            }
        } else {
            if (aggregate == "aminoAcid") {
                merged <- aggregate(count ~ aminoAcid, data = productive.seqs, FUN = sum)
                merged$frequencyCount <- merged$count/sum(merged$count) * 100
            }
            if (aggregate == "nucleotide") {
                count <- aggregate(count ~ nucleotide, data = productive.seqs, FUN = sum)
                productive.seqs$count <- NULL
                merged <- merge(productive.seqs, count, by = "nucleotide")
                merged$frequencyCount <- merged$count/sum(merged$count) * 100
            }
        }
        merged <- merged[order(merged$frequencyCount, decreasing = TRUE), intersect(names(sample), names(merged))]
        rownames(merged) <- NULL
        return(merged)
    }
} 