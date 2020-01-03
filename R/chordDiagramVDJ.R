#' Chord diagram of VJ or DJ gene associations
#' 
#' Creates a chord diagram showing VJ or DJ gene associations from one or more 
#' samples.
#' 
#' @param sample A data frame consisting of frequencies of antigen receptor 
#' sequences.  "vFamily", "jFamily", and if applicable, "dFamily" 
#' are a required columns.  Using output from the LymphoSeq function topSeqs is
#' recommended.
#' @param association A character vector of gene familes to associate.  Options 
#' include "VJ" or "DJ".
#' @param colors A character vector of 2 colors corresponding to the V/D and J 
#' gene colors respectively.
#' @details The size of the ribbons connecting VJ or DJ genes correspond to the 
#' number of samples or number of sequences that make up that recombination 
#' event.  The thicker the ribbon, the higher the frequency of the recombination.
#' @return Returns a chord diagram showing VJ or DJ gene associations from one or 
#' more samples.
#' @seealso \code{\link{topSeqs}}
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.nt <- productiveSeq(file.list = file.list, aggregate = "nucleotide")
#' 
#' top.seqs <- topSeqs(productive.seqs = productive.nt, top = 1)
#' 
#' chordDiagramVDJ(sample = top.seqs, association = "VJ", colors = c("red", "blue"))
#' 
#' # Remove "TCRB" from gene family name
#' top.seqs <- as.data.frame(apply(top.seqs, 2, function(x) gsub("TCRB", "", x)))
#' 
#' chordDiagramVDJ(sample = top.seqs, association = "VJ", colors = c("red", "blue"))
#' @export
#' @importFrom circlize colorRamp2 chordDiagram
chordDiagramVDJ <- function(sample, association = "VJ", colors = c("red", "blue")) {
    if (association == "VJ") {
        if (!all(c("vFamily", "jFamily") %in% colnames(sample))) {
            stop("The source data frame does not contain the required columns 'vFamily' and 'jFamily'.")
        }
        vj <- sample[, c("vFamily", "jFamily")]
        vj[is.na(vj) | vj == ""] <- "Unresolved"
        vj$vFamily <- as.character(vj$vFamily)
        vj$jFamily <- as.character(vj$jFamily)
        table <- table(vj$vFamily, vj$jFamily)
        matrix <- as.matrix(as.data.frame.matrix(table))
        ribbon.color <- circlize::colorRamp2(range(matrix), c("grey", "black"))
        circlize::chordDiagram(matrix, 
                               annotationTrack = c("grid", "name"), 
                               grid.col = c(rep(colors[1], dim(matrix)[1]), rep(colors[2], dim(matrix)[2])), 
                               col = ribbon.color)
    }
    if (association == "DJ") {
        if (!all(c("dFamily", "jFamily") %in% colnames(sample))) {
            stop("The source data frame does not contain the required columns 'dFamily' and 'jFamily'.")
        }
        dj <- sample[, c("dFamily", "jFamily")]
        dj[is.na(dj) | dj == ""] <- "Unresolved"
        dj$dFamily <- as.character(dj$dFamily)
        dj$jFamily <- as.character(dj$jFamily)
        table <- table(dj$dFamily, dj$jFamily)
        matrix <- as.matrix(as.data.frame.matrix(table))
        ribbon.color <- circlize::colorRamp2(range(matrix), c("grey", "black"))
        circlize::chordDiagram(matrix, 
                               annotationTrack = c("grid", "name"), 
                               grid.col = c(rep(colors[1], dim(matrix)[1]), rep(colors[2], dim(matrix)[2])), 
                               col = ribbon.color)
    }
}