#' Productive sequences
#' 
#' Remove unproductive CDR3 sequences from a list of data frames.
#' 
#' @param file.list A list of data frames consisting antigen receptor sequencing 
#' data imported by the LymphoSeq function readImmunoSeq. "aminoAcid", "count", and 
#' "frequencyCount" are required columns.
#' @param aggregate Indicates whether the values of "count", "frequencyCount", 
#' and "esimatedNumberGenomes" should be aggregated by amino acid or nucleotide 
#' sequence.  Acceptable values are "aminoAcid" or "nucleotide".  If "aminoAcid" 
#' is selected, then resulting data frame will have columns corresponding to 
#' aminoAcid, count, frequencyCount, and estimatedNumberGenomes (if this 
#' column is available).  If "nucleotide" is selected then all columns in the 
#' original list will be present in the outputted list.  The difference in 
#' output is due to the fact that the same amino acid CDR3 sequence may be 
#' encoded by multiple unique nucleotide sequences with differing V, D, and J 
#' genes.
#' @param prevalence A Boolean value indicating if a new column should be added 
#' to each of the data frames giving the prevalence of each CDR3 amino acid 
#' sequence in 55 healthy donor peripheral blood samples.  TRUE means the column 
#' is added and FALSE means it is not.  Values range from 0 to 100\% where 
#' 100\% means the sequence appeared in the blood of all 55 individuals.  The 
#' prevalenceTRB database is located in a separate package called LymphoSeqDB 
#' that should be loaded automatically.
#' @return Returns a list of data frames of productive amino acid sequences with
#' recomputed values for "count", "frequencyCount", and 
#' "esimatedNumberGenomes".  A productive sequences is defined as a sequences 
#' that is in frame and does not have an early stop codon.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.nt <- productiveSeq(file.list = file.list, 
#'    aggregate = "nucleotide", prevalence = FALSE)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, 
#'   aggregate = "aminoAcid", prevalence = TRUE)
#' @seealso Refer to the LymphoSeqDB package for details regarding the 
#' prevalenceTRB database.
#' @export
#' @importFrom plyr llply
#' @import LymphoSeqDB
productiveSeq <- function(file.list, aggregate = "aminoAcid", prevalence = FALSE) {
    productive.list <- plyr::llply(file.list, function(x) 
        productive(x, aggregate), .progress = "text")
    if (prevalence == TRUE) {
        if (aggregate != "aminoAcid") {
            stop("In order to add prevalence to your list of data frames, aggregate must be equal 'aminoAcid'.", call. = FALSE)
        }
        cat("Adding prevalence column...")
        i <- 1
        productive.list.prevalence <- list()
        for(i in 1:length(productive.list)) {
            file <- productive.list[[i]]
            file$prevalence <- LymphoSeqDB::prevalenceTRB[match(file$aminoAcid, LymphoSeqDB::prevalenceTRB$aminoAcid), "prevalence"]
            file$prevalence[is.na(file$prevalence)] <- 0
            productive.list.prevalence <- c(productive.list.prevalence, list(file))
        }
        names(productive.list.prevalence) <- names(productive.list)
        return(productive.list.prevalence)
    } else {
        return(productive.list)
    }
} 