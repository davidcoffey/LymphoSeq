#' Align mutliple sequences
#' 
#' Perform multiple sequence alignment using one of three methods and output results as a pdf file.
#' 
#' @param list A list of data frames consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param sample A character vector indicating the name of the sample in the productive
#' sequence list. 
#' @param type A character vector indicating whether "aminoAcid" or "nucleotide" sequences
#' should be aligned.  If "aminoAcid" is specified, then run productiveSeqs first.
#' @param method A character vector indicating the multiple sequence alignment method to 
#' be used.  Refer to the Bioconductor msa package for more details.  Options incude 
#' "ClustalW", "ClustalOmega", and "Muscle".
#' @return Saves a pdf file of the multiple sequence alignemnt to the working directory.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.nt <- productiveSeq(file.list = file.list, aggregate = "nucleotide")
#' 
#' alignSeq(list = productive.nt, sample = "TRB_Unsorted_1320", type = "nucleotide", 
#'          method = "ClustalW")
#' @export
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @import msa
alignSeq = function(list, sample, type = "nucleotide", method = "ClustalOmega") {
    file = list[[sample]]
    if(nrow(file) > 150){
        file = file[1:150,]
        message("Only the first 150 sequences will be aligned")
    }
    file[file == ""] <- "unresolved"
    if(type == "nucleotide"){
        file = file[nchar(file$nucleotide) > 9, ]
        string = Biostrings::DNAStringSet(file$nucleotide)
    }
    if(type == "aminoAcid"){
        file = file[nchar(file$aminoAcid) > 3, ]
        string = Biostrings::AAStringSet(file$aminoAcid)
    }
        names(string) = paste(file$vFamilyName, file$dFamilyName, file$jFamilyName, file$count)
        names(string) = gsub(names(string), pattern = "IGH|IGL|IGK|TCRB|TCRA", replacement = "")
        names(string) = gsub(names(string), pattern = "unresolved", replacement = "UNR")
        alignment = msa::msa(string, method = method)
        msa::msaPrettyPrint(alignment, output = "pdf", file = paste(names(list[sample]), ".pdf", sep = ""), 
                       showNumbering = "right", showNames = "left", showConsensus = "top",  
                       shadingMode = "similar", shadingColors = "blues",
                       paperWidth = ncol(alignment)*0.1 + 5, paperHeight = nrow(alignment)*0.1 + 5,
                       showLogo = "bottom", showLogoScale = "right", logoColors = "rasmol",
                       psFonts = TRUE, showLegend = TRUE, askForOverwrite = FALSE,
                       furtherCode=c("\\defconsensus{.}{lower}{upper}","\\showruler{1}{top}"))
        message(paste("Aligment files saved to", getwd()))
}