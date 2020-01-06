#' Read ImmunoSeq files
#' 
#' Imports tab-separated value (.tsv) files exported by the Adaptive 
#' Biotechnologies ImmunoSEQ analyzer and stores them as a list of data frames.
#' 
#' @param path Path to the directory containing tab-delimited files.  Only 
#' files with the extension .tsv are imported.  The names of the data frames are 
#' the same as names of the files.
#' @param columns Column names from the tab-delimited files that you desire to 
#' import, all others will be disregarded.  May use "all" to import all columns.  
#' A warning may be called if columns contain no data or must be coereced to a 
#' different class.  Usually this warning can be ignored.
#' @param recursive A Boolean value indicating whether tab-delimited files in 
#' subdirectories should be imported.  If TRUE, then all files in the parent as 
#' well as the subdirectory are imported.  If FALSE, then only files in the 
#' parent directory are imported.
#' @details May import tab-delimited files containing antigen receptor 
#' sequencing from other sources (e.g. miTCR or miXCR) as long as the column 
#' names are the same as used by ImmunoSEQ files.  Available column headings in 
#' ImmunoSEQ files are:  
#'     
#'     adaptiveV1:
#'     
#'     column list:
#'     nucleotide,aminoAcid,count (reads),frequencyCount (%),cdr3Length,vMaxResolved,vFamilyName,
#' vGeneName,vGeneAllele,vFamilyTies,vGeneNameTies,vGeneAlleleTies,dMaxResolved,dFamilyName,
#' dGeneName,dGeneAllele,dFamilyTies,dGeneNameTies,dGeneAlleleTies,jMaxResolved,jFamilyName,
#' jGeneName,jGeneAllele,jFamilyTies,jGeneNameTies,jGeneAlleleTies,vDeletion,n1Insertion,
#' d5Deletion,d3Deletion,n2Insertion,jDeletion,vIndex,n1Index,dIndex,n2Index,jIndex,
#' estimatedNumberGenomes,sequenceStatus,cloneResolved,vOrphon,dOrphon,jOrphon,vFunction,
#' dFunction,jFunction,fractionNucleated,vAlignLength,vAlignSubstitutionCount,
#' vAlignSubstitutionIndexes,vAlignSubstitutionGeneThreePrimeIndexes,vSeqWithMutations
#' 
#' minimum required columns:
#'     nucleotide,aminoAcid,count (reads),frequencyCount (%), vFamilyName, vGeneName,
#' dFamilyName, dGeneName, jGeneName, jFamilyName, estimatedNumberGenomes, sequenceStatus
#' 
#' adaptiveV2:
#'     
#'     column list:
#'     nucleotide,aminoAcid,count (templates/reads),frequencyCount (%),cdr3Length,vMaxResolved,
#' vFamilyName,vGeneName,vGeneAllele,vFamilyTies,vGeneNameTies,vGeneAlleleTies,dMaxResolved,
#' dFamilyName,dGeneName,dGeneAllele,dFamilyTies,dGeneNameTies,dGeneAlleleTies,jMaxResolved,
#' jFamilyName,jGeneName,jGeneAllele,jFamilyTies,jGeneNameTies,jGeneAlleleTies,vDeletion,
#' n1Insertion,d5Deletion,d3Deletion,n2Insertion,jDeletion,vIndex,n1Index,dIndex,n2Index,
#' jIndex,estimatedNumberGenomes,sequenceStatus,cloneResolved,vOrphon,dOrphon,jOrphon,vFunction,
#' dFunction,jFunction,fractionNucleated,vAlignLength,vAlignSubstitutionCount,vAlignSubstitutionIndexes,
#' vAlignSubstitutionGeneThreePrimeIndexes,vSeqWithMutations
#' 
#' minimum required columns:
#'     nucleotide,aminoAcid,count (templates/reads),frequencyCount (%), vFamilyName, vGeneName,
#' dFamilyName, dGeneName, jGeneName, jFamilyName, estimatedNumberGenomes, sequenceStatus
#' 
#' 
#' adaptiveV3:
#'     
#'     column list:
#'     nucleotide,aminoAcid,count,(templates),frequencyCount,(%),cdr3Length,vMaxResolved,vFamilyName,
#' vGeneName,vGeneAllele,vFamilyTies,vGeneNameTies,vGeneAlleleTies,dMaxResolved,dFamilyName,dGeneName,
#' dGeneAllele,dFamilyTies,dGeneNameTies,dGeneAlleleTies,jMaxResolved,jFamilyName,jGeneName,
#' jGeneAllele,jFamilyTies,jGeneNameTies,jGeneAlleleTies,vDeletion,n1Insertion,d5Deletion,d3Deletion,
#' n2Insertion,jDeletion,vIndex,n1Index,dIndex,n2Index,jIndex,estimatedNumberGenomes,sequenceStatus,
#' cloneResolved,vOrphon,dOrphon,jOrphon,vFunction,dFunction,jFunction,fractionNucleated
#' 
#' minimum required columns:
#'     nucleotide,aminoAcid,count (templates),frequencyCount (%), vFamilyName, vGeneName,
#' dFamilyName, dGeneName, jGeneName, jFamilyName, estimatedNumberGenomes, sequenceStatus
#' 
#' 
#' adaptiveV4:
#'     
#'     column list:
#'     nucleotide.CDR3.in.lowercase.,aminoAcid.CDR3.in.lowercase.,cloneCount,clonefrequency....,
#' CDR3Length,vGene,dGene,jGene,vDeletion,d5Deletion,d3Deletion,jDeletion,vdInsertion,
#' djInsertion,vjInsertion,fuction,CDR3.stripped.x.a
#' 
#' minimum required columns:
#'     nucleotide.CDR3.in.lowercase.,aminoAcid.CDR3.in.lowercase.,cloneCount,clonefrequency....,
#' vGene,dGene,jGene,fuction
#' 
#' bgiClone:
#'     
#'     column list:
#'     nucleotide(CDR3 in lowercase),aminoAcid(CDR3 in lowercase),cloneCount,clonefrequency (%),
#' CDR3Length,vGene,dGene,jGene,vDeletion,d5Deletion,d3Deletion,jDeletion,vdInsertion,
#' djInsertion,vjInsertion,fuction
#' 
#' minimum required columns:
#'     nucleotide(CDR3 in lowercase),aminoAcid(CDR3 in lowercase),cloneCount,clonefrequency (%),
#' vGene,dGene, fuction
#' @return Returns a list of data frames.  The names of each data frame are
#' assigned according to the original ImmunoSEQ file names.
#' @examples
#' file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")
#' 
#' file.list <- readImmunoSeq(path = file.path, 
#'                            columns = c("aminoAcid", "nucleotide", "count", 
#'                                      "count (templates)", "count (reads)", 
#'                                      "count (templates/reads)",
#'                                      "frequencyCount", "frequencyCount (%)", 
#'                                      "estimatedNumberGenomes"), 
#'                            recursive = FALSE)
#' @export
#' @importFrom data.table fread
#' @importFrom plyr llply 
#' @import tidyverse

getFileType <- function(clone_file) {
    adaptiveV1 <- c("nucleotide", "aminoAcid", "count (reads)", "frequencyCount (%)", "vGeneName", 
        "dGeneName", "jGeneName", "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
        "estimatedNumberGenomes")
        
    adaptiveV2 <- c("nucleotide", "aminoAcid", "count (templates/reads)", "frequencyCount (%)", "vGeneName", 
        "dGeneName", "jGeneName", "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
        "estimatedNumberGenomes")
    
    adaptiveV3 <- c("nucleotide", "aminoAcid", "count (templates)", "frequencyCount (%)", "vGeneName", 
        "dGeneName", "jGeneName", "vFamilyName", "dFamilyName", "jFamilyName", "sequenceStatus", 
        "estimatedNumberGenomes")
    
    adaptiveV4 <- c("nucleotide.CDR3.in.lowercase.", "aminoAcid.CDR3.in.lowercase.", "cloneCount", 
        "clonefrequency....", "vGene", "dGene", "jGene", "fuction")

    bgiClone <- c("nucleotide(CDR3 in lowercase)", "aminoAcid(CDR3 in lowercase)", "cloneCount",
        "clonefrequency (%)", "vGene", "dGene", "jGene", "fuction")

    columns <- invisible(colnames(readr::read_tsv(clone_file, n_max=1, col_types = cols())))
    if (all(adaptiveV1 %in% columns)) {
        file_type <- "adaptiveV1"
        header_list <- cols_only(nucleotide = "c", aminoAcid = "c", `count (reads)` = "i", 
                `frequencyCount (%)` = "d", vGeneName = "c", dGeneName = "c", 
                jGeneName = "c", vFamilyName = "c", dFamilyName = "c", 
                jFamilyName = "c", sequenceStatus = "c", estimatedNumberGenomes = "i")
    } else if (all(adaptiveV2 %in% columns)) {
        file_type <- "adaptiveV2"
        header_list <- cols_only(nucleotide = "c", aminoAcid = "c", `count (templates/reads)` = "i", 
                `frequencyCount (%)` = "d", vGeneName = "c", dGeneName = "c", 
                jGeneName = "c", vFamilyName = "c", dFamilyName = "c", 
                jFamilyName = "c", sequenceStatus = "c", estimatedNumberGenomes = "i")
    } else if (all(adaptiveV3 %in% columns)) {
        file_type <- "adaptiveV3"
        header_list <- cols_only(nucleotide = "c", aminoAcid = "c", `count (templates)` = "i", 
                `frequencyCount (%)` = "d", vGeneName = "c", dGeneName = "c", 
                jGeneName = "c", vFamilyName = "c", dFamilyName = "c", 
                jFamilyName = "c", sequenceStatus = "c", estimatedNumberGenomes = "i")
    } else if (all(adaptiveV4 %in% columns)) {
        file_type <- "adaptiveV4"
        header_list <- cols_only(`nucleotide.CDR3.in.lowercase.` = col_character(), 
                `aminoAcid.CDR3.in.lowercase.` = col_character(), cloneCount = col_integer(), 
                `clonefrequency....` = col_double(), vGene = col_character(), dGene = col_character(), 
                jGene = col_character(), fuction = col_character())
    } else if (all(bgiClone %in% columns)) {
        file_type <- "bgiClone"
        header_list <- cols_only(`nucleotide(CDR3 in lowercase)` = col_character(), 
                `aminoAcid(CDR3 in lowercase)` = col_character(), cloneCount = col_integer(), 
                `clonefrequency (%)` = col_double(), vGene = col_character(), 
                dGene = col_character(), jGene = col_character(), fuction = col_character())
    }
    ret_val <- list(file_type, header_list)
    return(ret_val)
}

getStandard <- function(file_type) {
    adaptiveV1 <- c(count = "count (reads)", frequencyCount = "frequencyCount (%)",
            `function` = "sequenceStatus")

    adaptiveV2 <- c(count = "count (templates/reads)", frequencyCount = "frequencyCount (%)",
            `function` = "sequenceStatus")

    adaptiveV3 <- c(count = "count (templates)", frequencyCount = "frequencyCount (%)",
            `function` = "sequenceStatus")

    adaptiveV4 <- c(nucleotide = "nucleotide.CDR3.in.lowercase.", aminoAcid = "aminoAcid.CDR3.in.lowercase.", 
            count = "cloneCount", frequencyCount = "clonefrequency....", vGeneName = "vGene", dGeneName = "dGene", 
            jGeneName = "jGene", `function` = "fuction")   
    
    bgiClone <- c(nucleotide = "nucleotide(CDR3 in lowercase)", aminoAcid = "aminoAcid(CDR3 in lowercase)", 
            count = "cloneCount", frequencyCount = "clonefrequency (%)", vGeneName = "vGene", dGeneName = "dGene", 
            jGeneName = "jGene", `function` = "fuction")

    type_hash <- list("adaptiveV1"=adaptiveV1, "adaptiveV2"=adaptiveV2, "adaptiveV3"=adaptiveV3, "adaptiveV4"=adaptiveV4,
            "bgiClone"=bgiClone)
    return(type_hash[[file_type]])
}

readFiles <- function(clone_file) {
    file_info <- getFileType(clone_file)
    file_type <- file_info[[1]]
    header_list <- file_info[[2]]
    col_std <- getStandard(file_type)
    col_old = header_list
    file_names <- tools::file_path_sans_ext(basename(clone_file))
    clone_frame <- readr::read_tsv(clone_file, col_types = col_old,
                                    na = c("", "NA", "Nan", "NaN"), trim_ws = TRUE,
                                    progress = show_progress())
    if ((file_type == "adaptiveV1") | (file_type == "adaptiveV2") | (file_type == "adaptiveV3")){
        clone_frame <- clone_frame %>% dplyr::rename(!!!col_std)
        clone_frame <- clone_frame %>% mutate(`function`= str_replace(`function`, "In", "in-frame"))
        clone_frame <- clone_frame %>% add_column(sample=file_names)
    } else if ((file_type == "adaptiveV4") | (file_type == "bgiClone")) {
        clone_frame <- clone_frame %>% dplyr::rename(!!!col_std)
        clone_frame <- clone_frame %>% extract(nucleotide, c("nucleotide"), regex = "([acgt]+)")
        clone_frame <- clone_frame %>% mutate(nucleotide  = toupper(nucleotide))
        clone_frame <- clone_frame %>% extract(aminoAcid, c("aminoAcid"), regex = "([acdefghiklmnpqrstvwy]+)")
        clone_frame <- clone_frame %>% mutate(aminoAcid = toupper(aminoAcid))
        clone_frame <- clone_frame %>% separate(vGeneName, c("vFamilyName"), sep = '-', remove = FALSE)
        clone_frame <- clone_frame %>% separate(dGeneName, c("dFamilyName"), sep = '-', remove = FALSE)
        clone_frame <- clone_frame %>% separate(jGeneName, c("jFamilyName"), sep = '-', remove = FALSE)
        clone_frame <- clone_frame %>% add_column(sample=file_names)

        clone_frame <- clone_frame %>% group_by(nucleotide) %>% 
            summarize(aminoAcid = first(aminoAcid), `count` = sum(`count`), vGeneName = first(vGeneName), `function` = first(`function`),
            jGeneName = first(jGeneName), dGeneName = first(dGeneName), vFamilyName = first(vFamilyName), dFamilyName = first(dFamilyName),
            jFamilyName = first(jFamilyName)) %>% mutate(frequencyCount = `count` / sum(`count`))
        clone_frame <- clone_frame %>% mutate(estimatedNumberGenomes=`count`)        
    }
    return(clone_frame)
}



readImmunoSeq <- function(path, recursive = FALSE) {
    file.paths1 <- list.files(path, full.names = TRUE, all.files = FALSE, 
                             recursive = recursive, pattern = ".tsv", 
                             include.dirs = FALSE)
    file.info <- file.info(file.paths1)
    file.paths2 <- rownames(file.info)[file.info$size > 0]
    if(!identical(file.paths1,file.paths2)){
        warning("One or more of the files you are trying to import has no sequences and will be ignored.", call. = FALSE)
    }
    file.list <- suppressWarnings(plyr::llply(file.paths2, readFiles, .progress = "text"))
    if(length(unique(plyr::llply(file.list, ncol))) > 1){
        warning("One or more of the files you are trying to import do not contain all the columns you specified.", call. = FALSE)
    }
    file.names <- gsub(".tsv", "", basename(file.paths2))
    names(file.list) <- file.names
    return(file.list)
} 