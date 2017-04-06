#' Differential abundance analysis
#' 
#' Use a Fisher exact test to calculate differential abdunance of each sequence in 
#' two samples and report the log2 transformed fold change, P value and adjusted P value.
#' 
#' @param list A list of data frames consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param sample1 A character vector indicating the name of the first sample in the list
#' to be compared.
#' @param sample2 A character vector indicating the name of the second sample in the list
#' to be compared.
#' @param type A character vector indicating whether "aminoAcid" or "nucleotide" sequences
#' should be used.  If "aminoAcid" is specified, then run productiveSeqs first.
#' @param q A numeric value between 0.0 and 1.0 indicating the threshold Holms adjusted 
#' P value (also knowns as the false discovery rate or q value) to subset the results with.  
#' Any sequences with a q value greater than this value will not be shown.
#' @param zero A numeric value to set all zero values to when calculating the log2
#' transformed fold change between samples 1 and 2.  This does not apply to the 
#' p and q value calculations.
#' @param abundance The input value for the Fisher exact test.  "estimatedNumberGenomes"
#' is recommend but "count" may also be used.
#' @param parallel A boolean indicating wheter parallel processing should be used or not.
#' @param cores An integer specifying the number of processor cores if parallel processing
#' is used.
#' @return Returns a data frame with columns corresponding to the frequency of the abudance
#' measure in samples 1 and 2, the P value, Q value (Holms adjusted P value, also knowns as 
#' the false discovery rate), and log2 transformed fold change.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' differentialAbundance(list = productive.aa, sample1 = "TRB_Unsorted_949", 
#'                       sample2 = "TRB_Unsorted_1320", type = "aminoAcid", q = 0.01, 
#'                       zero = 0.001, parallel = FALSE)
#' @export
#' @importFrom stats p.adjust
#' @importFrom stats fisher.test
#' @importFrom doMC registerDoMC
differentialAbundance <- function(sample1, sample2, list, abundance = "estimatedNumberGenomes", 
                                  type = "aminoAcid", q = 1, zero = 0.001, parallel = TRUE, cores = 4) {
    x <- list[[sample1]]
    y <- list[[sample2]]
    sequences <- union(x[, type], y[, type])
    fisherFunction <- function(sequences) {
        p.value <- as.numeric()
        sample1.freq <- as.numeric()
        sample2.freq <- as.numeric()
        if (any(x[, type] == sequences)) {
            in.x <- x[x[, type] == sequences, abundance]
            x.freq <- x[x[, type] == sequences, abundance]/sum(x[, abundance]) * 
                100
        } else {
            in.x <- 0
            x.freq <- 0
        }
        if (any(y[, type] == sequences)) {
            in.y <- y[y[, type] == sequences, abundance]
            y.freq <- y[y[, type] == sequences, abundance]/sum(y[, abundance]) * 
                100
        } else {
            in.y <- 0
            y.freq <- 0
        }
        not.x <- sum(x[, abundance]) - in.x
        not.y <- sum(y[, abundance]) - in.y
        matrix <- matrix(c(in.x, in.y, not.x, not.y), nrow = 2)
        fisher <- stats::fisher.test(matrix, workspace = 2e6)
        p.value <- c(p.value, fisher$p.value)
        sample1.freq <- c(sample1.freq, x.freq)
        sample2.freq <- c(sample2.freq, y.freq)
        results <- data.frame(aminoAcid = sequences, 
                              x = sample1.freq, 
                              y = sample2.freq, 
                              p = p.value)
        return(results)
    }
    doMC::registerDoMC(cores = cores)
    results <- plyr::ldply(sequences, fisherFunction, .parallel = parallel)
    results$q <- stats::p.adjust(results$p, method = "holm")
    x.not.zero <- results[, "x"]
    x.not.zero[x.not.zero == 0] <- zero
    y.not.zero <- results[, "y"]
    y.not.zero[y.not.zero == 0] <- zero
    results$l2fc <- log2(x.not.zero/y.not.zero)
    names(results)[2] <- sample1
    names(results)[3] <- sample2
    results <- results[order(results$q, decreasing = FALSE), ]
    results <- results[results$q <= q, ]
    rownames(results) <- NULL
    return(results)
}