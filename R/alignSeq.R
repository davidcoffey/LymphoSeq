#' Align mutliple sequences
#' 
#' Perform multiple sequence alignment using one of three methods and output results to the 
#' console or as a pdf file.  One may perform the alignment of all amino acid or nucleotide
#' sequences in a single sample.  Alternatively, one may search for a given sequence 
#' within a list of samples using an edit distance threshold.
#' 
#' @param list A list of data frames consisting of antigen receptor sequences 
#' imported by the LymphoSeq function readImmunoSeq.
#' @param sample A character vector indicating the name of the sample in the productive
#' sequence list. 
#' @param sequence A character vector of one ore more amino acid or nucleotide 
#' CDR3 sequences to search.
#' @param editDistance An integer giving the minimum edit distance that the 
#' sequence must be less than or equal to.  See details below.
#' @param type A character vector indicating whether "aminoAcid" or "nucleotide" sequences
#' should be aligned.  If "aminoAcid" is specified, then run productiveSeqs first.
#' @param method A character vector indicating the multiple sequence alignment method to 
#' be used.  Refer to the Bioconductor msa package for more details.  Options incude 
#' "ClustalW", "ClustalOmega", and "Muscle".
#' @param output A character vector indicating where the multiple sequence alignemnt should be
#' printed.  Options include "console" or "pdf".  If "pdf" is selected, the file is saved to
#' the working directory.  For "pdf" to work, Texshade must be installed.  Refer to the 
#' Bioconductor package msa installation instructions for more details.
#' @details Edit distance is a way of quantifying how dissimilar two sequences 
#' are to one another by counting the minimum number of operations required to 
#' transform one sequence into the other.  For example, an edit distance of 0 
#' means the sequences are identical and an edit distance of 1 indicates that 
#' the sequences different by a single amino acid or nucleotide.
#' 
#' @return Performs a multiple sequence alignemnt and prints to the console or saves a pdf to 
#' the working directory.
#' @seealso If having trouble saving pdf files, refer to Biconductor package msa for
#' installation instructions
#' \url{http://bioconductor.org/packages/release/bioc/vignettes/msa/inst/doc/msa.pdf}
#' @examples
#' file.path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path)
#' 
#' productive.nt <- productiveSeq(file.list = file.list, aggregate = "nucleotide")
#' 
#' alignSeq(list = productive.nt, sample = "IGH_MVQ92552A_BL", type = "nucleotide", 
#'          method = "ClustalW", output = "console")
#' @export
#' @importFrom Biostrings DNAStringSet
#' @importFrom Biostrings AAStringSet
#' @import msa
alignSeq = function(list, sample = NULL, sequence = NULL, editDistance = 15, output = "console", type = "nucleotide", method = "ClustalOmega") {
    if(!is.null(sequence) & is.null(sample)){
        file <- searchSeq(list = list, sequence = sequence, editDistance = editDistance, type = type, match = "partial")
        if(is.null(file)){
            stop("There are no sequences to be aligned", call. = FALSE)
        }
    }
    if(!is.null(sequence) & !is.null(sample)){
        file <- searchSeq(list = list[sample], sequence = sequence, editDistance = editDistance, type = type, match = "partial")
        if(is.null(file)){
            stop("There are no sequences to be aligned", call. = FALSE)
        }
    }
    if(is.null(sequence)){
        file <- list[[sample]] 
    }
    if(nrow(file) > 150){
        file <- file[1:150,]
        message("Only the first 150 sequences will be aligned")
    }
    file[file == ""] <- "unresolved"
    if(type == "nucleotide"){
        names(file)[which(names(file) == "foundSequnece")] <- "nucleotide"
        file <- file[nchar(file$nucleotide) > 3, ]
        string <- Biostrings::DNAStringSet(file$nucleotide)
    }
    if(type == "aminoAcid"){
        names(file)[which(names(file) == "foundSequnece")] <- "aminoAcid"
        file <- file[nchar(file$aminoAcid) > 3, ]
        string <- Biostrings::AAStringSet(file$aminoAcid)
    }
    if(nrow(file) < 3){
        stop("There are less than 3 sequences to be aligned", call. = FALSE)
    }
    if(!is.null(sequence)){
        names(string) <- paste(file$sample)
    } else {
        names(string) <- paste(file$vFamilyName, file$dFamilyName, file$jFamilyName, file$count)
    }
        names(string) <- gsub(names(string), pattern = "IGH|IGL|IGK|TCRB|TCRA", replacement = "")
        names(string) <- gsub(names(string), pattern = "unresolved", replacement = "UNR")
        alignment <- msa::msa(string, method = method)
        if(output == "console"){
            print(alignment, show = "complete")
        }
        if(output == "pdf"){
            if(!is.null(sample)){
                file.name <- paste(names(list[sample]), ".pdf", sep = "")
            } else {
                file.name <- "Sequence_aligment.pdf"
            }
            msa::msaPrettyPrint(alignment, output = "pdf", file = file.name, 
                                showNumbering = "right", showNames = "left", showConsensus = "top",  
                                shadingMode = "similar", shadingColors = "blues",
                                paperWidth = ncol(alignment)*0.1 + 5, paperHeight = nrow(alignment)*0.1 + 5,
                                showLogo = "bottom", showLogoScale = "right", logoColors = "rasmol",
                                psFonts = TRUE, showLegend = TRUE, askForOverwrite = FALSE,
                                furtherCode=c("\\defconsensus{.}{lower}{upper}","\\showruler{1}{top}"))
            message(paste(file.name, "saved to", getwd()))
        }
}
