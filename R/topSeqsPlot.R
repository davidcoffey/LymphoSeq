#' Cumulative frequency bar plot of top sequences
#' 
#' Create a cumulative frequency bar plot of a specified number of top 
#' sequences.
#' 
#' @param list A list data frames imported using the LymphoSeq function readImmunoSeq 
#' or productiveSeq.
#' @param top The number of top sequences to be colored in the bar plot.  All 
#' other, less frequent sequences are colored violet.
#' @return Returns a cumulative frequency bar plot of the top sequences.
#' @details The plot is made using the package ggplot2 and can be reformatted
#' using ggplot2 functions.  See examples below.
#' @seealso An excellent resource for examples on how to reformat a ggplot can 
#' be found in the R Graphics Cookbook online (\url{http://www.cookbook-r.com/Graphs/}).
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' topSeqsPlot(list = productive.aa, top = 10)
#' 
#' # Display the number of sequences at the top of bar plot and add a title
#' n <- as.character(lapply(file.list, nrow))
#' 
#' topSeqsPlot(list = productive.aa, top = 10) + 
#'    ggplot2::annotate("text", x = 1:length(file.list), y = 105, label = n, color = "black") +
#'    ggplot2::expand_limits(y = c(0, 110)) + ggplot2::ggtitle("Figure Title") + 
#'    ggplot2::scale_x_discrete(limits = names(file.list))
#' @export
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom plyr llply
topSeqsPlot <- function(list, top = 10) {
    if (any(top > lapply(list, nrow))) {
        stop(paste("The value for 'top' must be less than the smallest number of sequences in your data set (", min(unlist(lapply(list, nrow))), ")", sep = ""))
    }
    dominant <- plyr::llply(list, function(x) 
        x$frequencyCount[order(x$frequencyCount, decreasing = TRUE)][1:top])
    aminoAcid <- plyr::llply(list, function(x) 
        x$aminoAcid[order(x$frequencyCount, decreasing = TRUE)][1:top])
    subdominant <- plyr::llply(dominant, function(x) 100 - sum(x))

    dominant.df <- plyr::ldply(dominant, data.frame)
    aminoAcid.df <- plyr::ldply(aminoAcid, data.frame)
    subdominant.df <- plyr::ldply(subdominant, data.frame)
    
    dominant.df$aminoAcid <- aminoAcid.df$X..i..
    subdominant.df$aminoAcid <- rep("All other sequences", ncol(subdominant.df))
    topfreq <- rbind(dominant.df,subdominant.df)
    
    colnames(topfreq) <- c("Sample", "Frequency", "aminoAcid")
    topfreq$Sequence <- factor(paste("Sequence", c(rep(1:top, length(list)), 
                                                   rep(top + 1, length(list)))))
    topfreq$Sequence <- factor(topfreq$Sequence, 
                               levels = paste("Sequence", 1:(top + 1)))
    x.order <- names(subdominant[order(topfreq[topfreq$Sequence == paste("Sequence", top + 1), "Frequency"])])
    topfreq$CumulativeFrequency = topfreq$Frequency
    getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(11, "Spectral"))
    ggplot(topfreq, aes_string(x = "Sample", y = "CumulativeFrequency", fill = "Sequence", label = "Frequency", text = "aminoAcid")) + 
        geom_bar(stat = "identity") + 
        scale_x_discrete(limits = x.order) + 
        scale_fill_manual(values = getPalette(top + 1)) + 
        theme_classic() + 
        scale_y_continuous(expand = c(0, 0)) + 
        theme(legend.position = "none") + 
        theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1, size = 10), 
              axis.text.y = element_text(size = 10)) + labs(x = "", y = "Frequency (%)")
} 
