#' Pairwise comparison plot
#' 
#' Creates a heat map from a similarity or Bhattacharyya matrix.
#' 
#' @param matrix A similarity or Bhattacharyya matrix produced by the LymphoSeq 
#' functions similarityMatrix or bhattacharyyaMatrix.
#' @return A pairwise comparison heat map.
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can 
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' The functions to create the similarity or Bhattacharyya matrix can be found 
#' here: \code{\link{similarityMatrix}} and \code{\link{bhattacharyyaMatrix}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' similarity.matrix <- similarityMatrix(productive.seqs = productive.aa)
#' 
#' pairwisePlot(matrix = similarity.matrix)
#' 
#' bhattacharyya.matrix <- bhattacharyyaMatrix(productive.seqs = productive.aa)
#' 
#' pairwisePlot(matrix = bhattacharyya.matrix)
#' 
#' # Change plot color, title legend, and add title
#' pairwisePlot(matrix = similarity.matrix) + 
#'    ggplot2::scale_fill_gradient(low = "#deebf7", high = "#3182bd") + 
#'    ggplot2::labs(fill = "Similarity score") + ggplot2::ggtitle("Figure Title")
#' @export
#' @import ggplot2
#' @importFrom reshape melt.data.frame
#' @importFrom stats na.omit
pairwisePlot <- function(matrix) {
    i <- 1
    l <- length(matrix)
    p <- l - 1
    for (i in 1:p) {
        j <- i + 1
        matrix[i, j:l] <- NA
    }
    matrix$Samples <- rownames(matrix)
    melt <- reshape::melt.data.frame(matrix, id.vars = "Samples")
    melt <- na.omit(melt)
    melt$variable <- factor(melt$variable, levels = rownames(matrix))
    melt$Samples <- factor(melt$Samples, levels = rev(rownames(matrix)))
    names(melt) = c("Sample.x", "Sample.y", "Score")
    ggplot(data = melt, aes_string(x = "Sample.x", y = "Sample.y", fill = "Score")) + 
        geom_tile() + 
        scale_fill_gradient(low = "#fee8c8", high = "#e34a33") + 
        theme_classic() + 
        labs(x = "", y = "", fill = "") + 
        theme(axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5))
}