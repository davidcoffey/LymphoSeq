#' Common sequences Venn diagram
#' 
#' Creates a Venn diagram comparing the number of common sequences in two or 
#' three samples.
#' 
#' @param samples A character vector of two or three names of samples in 
#' productive.seqs to compare.
#' @param productive.seqs A list of productive amino acid sequences generated
#' by the LymphoSeq function productiveSeq.
#' @return Returns a a Venn diagram of the number of common sequences between
#' two or three samples.
#' @seealso \code{\link{commonSeqs}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.aa <- productiveSeq(file.list = file.list, aggregate = "aminoAcid")
#' 
#' # Plot a triple Venn diagram
#' commonSeqsVenn(samples = c("TCRB_Day0_Unsorted", 
#'    "TCRB_Day32_Unsorted", "TCRB_Day83_Unsorted"), 
#'    productive.seqs = productive.aa)
#' 
#' # Plot a double Venn diagram
#' commonSeqsVenn(samples = c("TCRB_Day0_Unsorted", 
#'    "TCRB_Day32_Unsorted"), productive.seqs = productive.aa)
#' 
#' # Save Venn diagram as a .png file to working directory
#' png(filename = "Venn diagram.png", res = 300, units = "in", height = 5, width = 5)
#' 
#' commonSeqsVenn(samples = c("TCRB_Day0_Unsorted", "TCRB_Day32_Unsorted"), 
#'    productive.seqs = productive.aa)
#' 
#' dev.off()
#' @export
#' @importFrom VennDiagram draw.pairwise.venn draw.triple.venn
#' @importFrom grid grid.newpage grid.draw
commonSeqsVenn <- function(samples, productive.seqs) {
    if (length(samples) > 3 | length(samples) < 2) {
        stop("Please enter 2 or 3 samples.")
    }
    if (length(samples) == 2) {
        a <- productive.seqs[[samples[1]]]
        b <- productive.seqs[[samples[2]]]
        grid::grid.newpage()
        venn <- VennDiagram::draw.pairwise.venn(area1 = length(a$aminoAcid), 
                                                area2 = length(b$aminoAcid), 
                                                cross.area = length(intersect(a$aminoAcid, b$aminoAcid)), 
                                                category = c(samples[1], samples[2]), 
                                                cat.fontfamily = rep("sans", 2), 
                                                fontfamily = rep("sans", 3), 
                                                fill = c("#3288bd", "#d53e4f"), 
                                                cat.pos = c(0, 0),
                                                cat.dist = rep(0.025, 2),
                                                cex = 1, 
                                                cat.cex = 0.7,
                                                lwd = rep(2, 2))
        grid::grid.draw(venn)
    }
    if (length(samples) == 3) {
        a <- productive.seqs[[samples[1]]]
        b <- productive.seqs[[samples[2]]]
        c <- productive.seqs[[samples[3]]]
        grid::grid.newpage()
        venn <- VennDiagram::draw.triple.venn(area1 = length(a$aminoAcid), 
                                              area2 = length(b$aminoAcid), 
                                              area3 = length(c$aminoAcid), 
                                              n12 = length(intersect(a$aminoAcid, b$aminoAcid)), 
                                              n23 = length(intersect(b$aminoAcid, c$aminoAcid)), 
                                              n13 = length(intersect(a$aminoAcid, c$aminoAcid)), 
                                              n123 = length(Reduce(intersect, list(a$aminoAcid, b$aminoAcid, c$aminoAcid))), 
                                              category = c(samples[1], samples[2], samples[3]), 
                                              cat.fontfamily = rep("sans", 3), 
                                              fontfamily = rep("sans", 7), 
                                              fill = c("#3288bd", "#abdda4", "#d53e4f"), 
                                              cat.pos = c(0, 0, 180), 
                                              cat.dist = rep(0.025, 3),
                                              cex = 1, 
                                              cat.cex = 0.7,
                                              lwd = rep(2, 3))
        grid::grid.draw(venn)
    }
} 