#' Unique sequences
#' 
#' Aggregates all productive sequences within a list of data frames by count.
#' 
#' @param productive.aa A list data frames of of productive amino acid sequences 
#' imported using the function LymphoSeq function productiveSeq where the 
#' aggregate parameter was set to "aminoAcid". 
#' @return A data frame of unique amino acid sequences from the list of 
#' data frames aggregated by count.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' unique.seqs <- uniqueSeqs(productive.aa = productive.aa)
#' @export
#' @importFrom plyr llply ldply
#' @importFrom stats aggregate
uniqueSeqs <- function(productive.aa) {
    if(any(unlist(lapply(productive.aa, function(x) 
        x[, "aminoAcid"] == "" |
        grepl("\\*", x[, "aminoAcid"]) | 
        duplicated(x[, "aminoAcid"]))))){
        stop("Your list contains unproductive sequences or has not been aggreated for productive amino acid sequences.  Remove unproductive sequences first using the function productiveSeq with the aggregate parameter set to 'aminoAcid'.", call. = FALSE)
    }
    unique.seq <- plyr::ldply(productive.aa, function(x) x[, c("aminoAcid", "count")])
    unique.seq$.id <- NULL
    unique.seq <- aggregate(data = unique.seq, count ~ aminoAcid, sum)
    unique.seq <- unique.seq[order(unique.seq$count, decreasing = TRUE),]
    return(unique.seq)
}
