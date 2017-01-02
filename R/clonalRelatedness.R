#' Clonal relatedness
#' 
#' Calculates the clonal relatedness for each sample in a list of data frames.
#' 
#' @param list A list data frames of unproductive or productive nucleotide 
#' sequences or productive nucleotide sequences.  Nucleotide and count are 
#' required columns.
#' @param editDistance An integer giving the minimum edit distance that the 
#' sequence must be less than or equal to. See details below.
#' @details  Clonal relatedness is the proportion of nucleotide sequences that
#' are related by a defined edit distance threshold.  The value ranges from 0 to
#' 1 where 0 indicates no sequences are related and 1 indicates all sequences 
#' are related.
#' 
#' Edit distance is a way of quantifying how dissimilar two sequences 
#' are to one another by counting the minimum number of operations required to 
#' transform one sequence into the other. For example, an edit distance of 0 
#' means the sequences are identical and an edit distance of 1 indicates that 
#' @return Returns a data frame with the calculated clonal relatedness for each sample.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' clonal.relatedness <- clonalRelatedness(list = file.list, editDistance = 10)
#' @export
#' @importFrom stringdist stringdist
clonalRelatedness <- function(list, editDistance = 10){
    relatedness <- as.numeric()
    clonality <- as.numeric()
    cloneResolved <- as.character()
    i <- 2
    for(i in 1:length(list)){
        file <- list[[i]]
        top <- file[order(file$count, decreasing = TRUE), "nucleotide"][1]
        result1 <- length(which(
            stringdist::stringdist(top, file$nucleotide, 
                                   method = "lv") < editDistance))/nrow(file)
        relatedness <- c(relatedness, result1)
    }
    df <- data.frame(samples = names(list), clonalRelatedness = relatedness)
    return(df)
}
