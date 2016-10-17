#' Top frequencies
#' 
#' Creates a data frame of the top productive amino acid sequences that have a 
#' specified minimum frequency threshold and reports the number of samples that 
#' the sequence appears in along with the minimum, maximum, and mean frequency 
#' across all samples.  For T cell receptor beta sequences, the \% prevalence 
#' and antigen specificity of that sequence is also provided.
#' 
#' @param productive.aa A list data frames of of productive amino acid sequences 
#' imported using the function LymphoSeq function productiveSeq where the 
#' aggregate parameter was set to "aminoAcid".  
#' @param percent The minimum \% frequency that the sequence appears in any of 
#' the listed samples.
#' @return A data frame of amino acid sequences and the number of samples that 
#' the sequence appears in along with the minimum, maximum, and mean frequency 
#' across all samples.  
#' For T cell receptor beta sequences, additionally reported is the 
#' \% prevalence that the sequence appears in 55 healthy donor blood samples.  
#' Also provided is the antigen specificity of that sequence if known by 
#' comparing it to a database of previously reported sequences in the 
#' literature.  The prevalenceTRB and publishedTRB databases are located in a 
#' separate package called LymphoSeqDB that should be loaded automatically.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' top.freq <- topFreq(productive.aa = productive.aa, percent = 0.1)
#' @seealso Refer to the LymphoSeqDB package for details regarding the 
#' prevalenceTRB and publishedTRB database.
#' @export
#' @importFrom dplyr group_by summarise
#' @import LymphoSeqDB
#' @import LymphoSeqDB
topFreq <- function(productive.aa, percent = 0.1) {
    if(any(unlist(lapply(productive.aa, function(x) 
        x[, "aminoAcid"] == "" | grepl("\\*", x[, "aminoAcid"]) | 
        duplicated(x[, "aminoAcid"]))))){
        stop("Your list contains unproductive sequences or has not been aggreated for productive amino acid sequences.  Remove unproductive sequences first using the function productiveSeq with the aggregate parameter set to 'aminoAcid'.", call. = FALSE)
    }
    frequencyCount <- NULL
    aminoAcid <- NULL
    top.seqs <- plyr::llply(productive.aa, function(x) 
        subset(x, frequencyCount >= percent, 
               select = aminoAcid))
    top.seqs <- Reduce(rbind, top.seqs)
    top.merged <- Reduce(rbind, productive.aa)
    top.merged <- top.merged[top.merged$aminoAcid %in% unique(top.seqs$aminoAcid), ]
    top.grouped <- dplyr::group_by(top.merged, aminoAcid)
    top.table <- data.frame(dplyr::summarise(top.grouped, 
                                             min(frequencyCount), 
                                             max(frequencyCount), 
                                             mean(frequencyCount), 
                                             length(aminoAcid)))
    colnames(top.table) <- c("aminoAcid", "minFrequency", "maxFrequency", 
                             "meanFrequency", "numberSamples")
    top.table <- top.table[order(top.table$numberSamples, 
                                 top.table$meanFrequency, 
                                 decreasing = TRUE), ]
    top.table$prevalence <- LymphoSeqDB::prevalenceTRB[match(top.table$aminoAcid, 
                                                             LymphoSeqDB::prevalenceTRB$aminoAcid), 
                                                       "prevalence"]
    top.table$prevalence[is.na(top.table$prevalence)] <- 0
    top.table$antigen <- as.character(LymphoSeqDB::publishedTRB[match(top.table$aminoAcid, 
                                                                      LymphoSeqDB::publishedTRB$aminoAcid), 
                                                                "antigen"])
    top.table <- replace(top.table, is.na(top.table), "")
    return(top.table)
}