#' Clone tracking plot
#' 
#' Creates line plot tracking amino acid frequencies across multiple samples
#' 
#' @param sequence.matrix A sequence matrix generated from the LymphoSeq 
#' function seqMatrix.
#' @param map An optional character vector of one or more sample names contained 
#' in the productive.aa list.  If the same sequence appears in multiple mapped 
#' samples, then it will be assigned the label of the first listed sample only.
#' @param productive.aa A list of data frames of productive amino acid sequences 
#' containing the samples to be mapped.  This parameter is only required if 
#' sequence mapping is performed.
#' @param label An optional character vector of one or more labels used to 
#' annotate the mapped sequences.  The order of the labels must match the order 
#' of the samples listed in map.
#' @param track An optional character vector of one or more amino acid sequences 
#' to track.
#' @param unassigned A Boolean value indicating whether or not to draw the lines 
#' of sequences not being mapped or tracked.  If TRUE then the unassigned 
#' sequences are drawn.  If FALSE, then the unassigned sequences are not drawn.
#' @return Returns a line plot showing the amino acid frequencies across 
#' multiple samples in the sequence matrix where each line represents one 
#' unique sequence.
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
#' top.freq <- topFreq(productive.aa = productive.aa, percent = 0.1)
#' 
#' sequence.matrix <- seqMatrix(productive.aa = productive.aa, sequences = top.freq$aminoAcid)
#' 
#' # Track clones without mapping or tracking specific sequences
#' cloneTrack(sequence.matrix = sequence.matrix)
#' 
#' # Track top 20 clones mapping to the CD4 and CD8 samples
#' cloneTrack(sequence.matrix = sequence.matrix, productive.aa = productive.aa, 
#'    map = c("TCRB_Day949_CD4", "TCRB_Day949_CD8"), label = c("CD4", "CD8"), 
#'    track = top.freq$aminoAcid[1:20], unassigned = TRUE) 
#' 
#' # Track the top 10 clones from top.freq
#' cloneTrack(sequence.matrix = sequence.matrix, productive.aa = productive.aa, 
#'    track = top.freq$aminoAcid[1:10], unassigned = FALSE) 
#' 
#' # Track clones mapping to the CD4 and CD8 samples while ignoring all others
#' cloneTrack(sequence.matrix = sequence.matrix, productive.aa = productive.aa, 
#'    map = c("TCRB_Day949_CD4", "TCRB_Day949_CD8"), label = c("CD4", "CD8"), 
#'    unassigned = FALSE) 
#' 
#' # Track clones mapping to the CD4 and CD8 samples and track 2 specific sequences
#' cloneTrack(sequence.matrix = sequence.matrix, productive.aa = productive.aa, 
#'    map = c("TCRB_Day949_CD4", "TCRB_Day949_CD8"), label = c("CD4", "CD8"), 
#'    track = c("CASSPPTGERDTQYF", "CASSQDRTGQYGYTF"), unassigned = FALSE)
#' 
#' # Reorder the x axis, change the axis labels, convert to log scale, and add title
#' x.limits <- c("TCRB_Day0_Unsorted", "TCRB_Day32_Unsorted", 
#'    "TCRB_Day83_Unsorted", "TCRB_Day949_Unsorted", "TCRB_Day1320_Unsorted")
#' 
#' sequence.matrix <- sequence.matrix[ ,c("aminoAcid", x.limits)]
#'    
#' clone.track <- cloneTrack(sequence.matrix = sequence.matrix, 
#'    productive.aa = productive.aa, track = top.freq$aminoAcid[1:10], unassigned = FALSE) 
#' 
#' x.labels <- c("Day 0", "Day 32", "Day 83", "Day 949", "Day 1320")
#' 
#' clone.track + 
#'    ggplot2::scale_x_discrete(expand = c(0,0), labels = x.labels) + 
#'    ggplot2::scale_y_log10() + ggplot2::annotation_logticks(sides = "l") + 
#'    ggplot2::ggtitle("Figure Title")
#' @export
#' @import ggplot2
#' @importFrom RColorBrewer brewer.pal
#' @importFrom reshape melt.data.frame
cloneTrack <- function(sequence.matrix, map = "none", productive.aa, 
                       label = "none", track = "none", unassigned = TRUE) {
    if (any(map != "none")) {
        sequence.matrix$numberSamples <- NULL
        sequence.melt <- reshape::melt.data.frame(sequence.matrix, 
                                                  id.vars = "aminoAcid")
        sequence.melt$map <- rep("Unassigned", length(sequence.matrix$aminoAcid))
        if (any(track != "none")) {
            i <- 1
            for (i in 1:length(track)) {
                sequence.melt$map <- ifelse(sequence.melt$aminoAcid == track[i], 
                                            track[i], sequence.melt$map)
            }
        }
        j <- 1
        for (j in 1:length(map)) {
            file <- productive.aa[[map[j]]]
            aminoAcid <- as.character(file$aminoAcid)
            sequence.melt$map <- ifelse(sequence.melt$aminoAcid %in% aminoAcid & 
                                            sequence.melt$map == "Unassigned", 
                                        label[j], sequence.melt$map)
        }
        if (unassigned == FALSE) {
            sequence.melt <- sequence.melt[sequence.melt$map != "Unassigned", ]
        }
        getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
        plot <- ggplot2::ggplot(sequence.melt, 
                                aes_string(x = "variable", 
                                           y = "value", 
                                           group = "aminoAcid", 
                                           color = "map")) + 
            geom_line() + 
            theme_minimal() + 
            scale_x_discrete(expand = c(0, 0)) + 
            scale_y_continuous(expand = c(0, 0)) + 
            labs(x = "", y = "Frequency (%)", color = "") + 
            scale_colour_manual(values = c(getPalette(length(label) + 
                                                          length(track) + 1))) + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        return(plot)
    } else {
        sequence.matrix$numberSamples <- NULL
        sequence.melt <- reshape::melt.data.frame(sequence.matrix, id.vars = "aminoAcid")
        if (any(track != "none")) {
            sequence.melt$map <- rep("Unassigned", length(sequence.matrix$aminoAcid))
            k <- 1
            for (k in 1:length(track)) {
                sequence.melt$map <- ifelse(sequence.melt$aminoAcid == track[k], 
                                            track[k], sequence.melt$map)
            }
            if (unassigned == FALSE) {
                sequence.melt <- sequence.melt[sequence.melt$map != "Unassigned", 
                                               ]
            }
            getPalette <- grDevices::colorRampPalette(RColorBrewer::brewer.pal(9, "Set1"))
            plot <- ggplot2::ggplot(sequence.melt, aes_string(x = "variable", 
                                                              y = "value", 
                                                              group = "aminoAcid", 
                                                              color = "map")) + 
                geom_line() + theme_minimal() + 
                scale_x_discrete(expand = c(0,0)) + 
                scale_y_continuous(expand = c(0, 0)) + 
                labs(x = "", y = "Frequency (%)", color = "") + 
                scale_colour_manual(values = getPalette(length(track) + 1)) + 
                theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
            return(plot)
        }
        plot <- ggplot2::ggplot(sequence.melt, aes_string(x = "variable", 
                                                          y = "value", 
                                                          group = "aminoAcid")) + 
            geom_line() + 
            theme_minimal() + 
            scale_x_discrete(expand = c(0, 0)) + 
            scale_y_continuous(expand = c(0, 0)) + 
            labs(x = "", y = "Frequency (%)") + 
            theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
        return(plot)
    }
} 