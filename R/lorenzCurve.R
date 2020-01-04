#' Lorenz curve
#' 
#' Plots a Lorenz curve derived from the frequency of the amino acid sequences.
#' 
#' @param samples A character vector of sample names in list.
#' @param list A list data frames generated using the LymphoSeq function readImmunoSeq 
#' or productiveSeq.  "frequencyCount" is a required column.
#' @return Returns a Lorenz curve.
#' @details The Gini coefficient is an alternative metric used to calculate 
#' repertoire diversity and is derived from the Lorenz curve.  The Lorenz curve 
#' is drawn such that x-axis represents the cumulative percentage of unique 
#' sequences and the y-axis represents the cumulative percentage of reads.  A 
#' line passing through the origin with a slope of 1 reflects equal frequencies 
#' of all sequences.  The Gini coefficient is the ratio of the area between the 
#' line of equality and the observed Lorenz curve over the total area under the 
#' line of equality.
#' 
#' The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can 
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' lorenzCurve(samples = names(file.list), list = file.list)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' lorenzCurve(samples = names(productive.aa), list = productive.aa)
#' 
#' # Change the legend labels, line colors, and add a title
#' samples <- c("TRB_Unsorted_0", "TRB_Unsorted_32", 
#'    "TRB_Unsorted_83", "TRB_Unsorted_949", "TRB_Unsorted_1320")
#' 
#' lorenz.curve <- lorenzCurve(samples = samples, list = productive.aa)
#' 
#' labels <- c("Day 0", "Day 32", "Day 83", "Day 949", "Day 1320")
#' 
#' colors <- c("navyblue", "red", "darkgreen", "orange", "purple")
#' 
#' lorenz.curve + ggplot2::scale_color_manual(name = "Samples", breaks = samples, 
#'    labels = labels, values = colors) + ggplot2::ggtitle("Figure Title")
#' @export
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom ineq Lc
lorenzCurve <- function(samples, list) {
    lorenz <- data.frame()
    i <- 1
    for (i in 1:length(samples)) {
        sample <- samples[i]
        file <- list[[sample]]
        lc <- ineq::Lc(file$frequencyCount)
        lcdf <- data.frame(L = lc$L, p = lc$p)
        lcdf$sample <- rep(sample, nrow(lcdf))
        lorenz <- rbind(lcdf, lorenz)
    }
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
    plot <- ggplot2::ggplot(lorenz, aes_string(x = "p", y = "L", color = "sample")) + 
        geom_line(size = 1) + 
        theme_minimal() + 
        scale_color_manual(values = getPalette(length(samples) + 1)) + 
        scale_y_continuous(expand = c(0, 0)) + 
        scale_x_continuous(expand = c(0,0)) + 
        geom_abline(intercept = 0, slope = 1, color = "grey", linetype = 2) + 
        coord_fixed() + 
        labs(x = "Cumulative percentage of unique sequences", 
             y = "Cumulative percentage of reads", color = "")
    return(plot)
}
