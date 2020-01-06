#' Chord diagram of VJ or DJ gene associations
#' 
#' Creates a chord diagram showing VJ or DJ gene associations from one or more 
#' samples.
#' 
#' @param sample A data frame consisting of frequencies of antigen receptor 
#' sequences.  "vFamilyName", "jFamilyName", and if applicable, "dFamilyName" 
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
        if (!all(c("vFamilyName", "jFamilyName") %in% colnames(sample))) {
            stop("The source data frame does not contain the required columns 'vFamilyName' and 'jFamilyName'.")
        }
        vj <- sample[, c("vFamilyName", "jFamilyName")]
        vj[is.na(vj) | vj == ""] <- "Unresolved"
        vj$vFamilyName <- as.character(vj$vFamilyName)
        vj$jFamilyName <- as.character(vj$jFamilyName)
        table <- table(vj$vFamilyName, vj$jFamilyName)
        matrix <- as.matrix(as.data.frame.matrix(table))
        ribbon.color <- circlize::colorRamp2(range(matrix), c("grey", "black"))
        circlize::chordDiagram(matrix, 
                               annotationTrack = c("grid", "name"), 
                               grid.col = c(rep(colors[1], dim(matrix)[1]), rep(colors[2], dim(matrix)[2])), 
                               col = ribbon.color)
    }
    if (association == "DJ") {
        if (!all(c("dFamilyName", "jFamilyName") %in% colnames(sample))) {
            stop("The source data frame does not contain the required columns 'dFamilyName' and 'jFamilyName'.")
        }
        dj <- sample[, c("dFamilyName", "jFamilyName")]
        dj[is.na(dj) | dj == ""] <- "Unresolved"
        dj$dFamilyName <- as.character(dj$dFamilyName)
        dj$jFamilyName <- as.character(dj$jFamilyName)
        table <- table(dj$dFamilyName, dj$jFamilyName)
        matrix <- as.matrix(as.data.frame.matrix(table))
        ribbon.color <- circlize::colorRamp2(range(matrix), c("grey", "black"))
        circlize::chordDiagram(matrix, 
                               annotationTrack = c("grid", "name"), 
                               grid.col = c(rep(colors[1], dim(matrix)[1]), rep(colors[2], dim(matrix)[2])), 
                               col = ribbon.color)
    }
}