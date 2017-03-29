#' Merge files
#' 
#' Merges two or more sample data frames into a single data frame and aggregates 
#' count, frequencyCount, and estimatedNumberGenomes.
#' 
#' @param samples A character vector of the names of the samples that are to be 
#' merged in file.list.
#' @param file.list A list of data frames imported using the LymphoSeq function
#' readImmunoSeq that contain the sample names that are to be merged.  "aminoAcid", 
#' "nucleotide", "count" and "frequencyCount" are required columns.
#' @return Returns a data frame of the merged samples.  The values of count, 
#' frequencyCount, and estimatedNumberGenomes are aggregated.  That is, the sum 
#' of count and estimatedNumberGenomes columns of the merged data frame should 
#' equal the sum of the columns from the unmerged samples.  Likewise, the 
#' frequencyCount of the merged data frame should add up to 100\%.  All other 
#' columns in the unmerged data frames are included in the merge data frame.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' TCRB_Day949_Merged <- mergeFiles(samples = c("TRB_CD4_949", 
#'    "TRB_CD8_949"), file.list)
#' 
#' # To combine the merged data frames with file.list
#' file.list <- c(list("TCRB_Day949_Merged" = TCRB_Day949_Merged), file.list)
#' @export
#' @importFrom stats aggregate
mergeFiles <- function(samples, file.list) {
    merged <- data.frame()
    i <- 1
    for (i in 1:length(samples)) {
        file <- file.list[[as.character(samples[i])]]
        merged <- rbind(merged, file)
    }
    count <- aggregate(count ~ nucleotide, data = merged, FUN = sum)
    if (any(grepl("estimatedNumberGenomes", colnames(merged)))) {
        merged$estimatedNumberGenomes[is.na(merged$estimatedNumberGenomes)] <- 0
        estimatedNumberGenomes <- aggregate(estimatedNumberGenomes ~ nucleotide, 
                                            data = merged, FUN = sum)
    }
    merged$count <- NULL
    merged$estimatedNumberGenomes <- NULL
    merged <- merged[!duplicated(merged$nucleotide), ]
    merged <- Reduce(function(x, y) 
        merge(x, y, by = "nucleotide"), list(count, estimatedNumberGenomes, merged))
    merged$frequencyCount <- (merged$count/sum(merged$count)) * 100
    return(merged)
} 