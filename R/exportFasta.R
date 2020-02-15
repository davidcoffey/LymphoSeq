#' Export sequences in fasta format
#' 
#' Export nucleotide or amino acid sequences in fasta format.
#' 
#' @param list A list of data frames consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param type A character vector indicating whether "aminoAcid" or "nucleotide" sequences
#' should be exported.  If "aminoAcid" is specified, then run productiveSeqs first.
#' @param names A character vector of one or more column names to name the sequences.
#' If "rank" is specified, then the rank order of the sequences by frequency is used.
#' @return Exports fasta files to the working directory.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' exportFasta(list = file.list, type = "nucleotide", names = c("rank", "aminoAcid", "count"))
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' exportFasta(list = productive.aa, type = "aminoAcid", names = "frequencyCount")
#' @export
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @importFrom Biostrings writeXStringSet
exportFasta <- function(list, type = "nucleotide", names = c("rank", "aminoAcid", "count")) {
    if (type == "nucleotide") {
        i <- 1
        for (i in 1:length(list)) {
            file <- list[[i]]
            sequences <- file$nucleotide
            fasta <- Biostrings::DNAStringSet(sequences)
            if (any(names %in% "rank")) {
                file$rank <- rownames(file)
            }
            names(fasta) <- do.call(paste, c(file[names], sep = "_"))
            Biostrings::writeXStringSet(fasta, paste(names(list[i]), ".fasta", sep = ""))
        }
    }
    if (type == "aminoAcid") {
        i <- 1
        for (i in 1:length(list)) {
            file <- list[[i]]
            sequences <- file$aminoAcid
            if (any(sequences == "") | any(grepl("\\*", sequences))) {
                stop("Your list contains unproductive sequences.  Remove unproductive sequences first using the function productiveSeq.", call. = FALSE)
            }
            fasta <- Biostrings::AAStringSet(sequences)
            if (any(names %in% "rank")) {
                file$rank <- rownames(file)
            }
            names(fasta) <- do.call(paste, c(file[names], sep = "_"))
            Biostrings::writeXStringSet(fasta, paste(names(list[i]), ".fasta", sep = ""))
        }
    }
    message(paste("Fasta files saved to", getwd()))
}
