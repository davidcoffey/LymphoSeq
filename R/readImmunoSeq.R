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
#' "nucleotide", "aminoAcid", "count", "count (templates)", "count (reads)", 
#' "count (templates/reads)", "frequencyCount", "frequencyCount (\%)", "cdr3Length", 
#' "vMaxResolved", "vFamilyName", "vGeneName", "vGeneAllele", "vFamilyTies", 
#' "vGeneNameTies", "vGeneAlleleTies", "dMaxResolved", "dFamilyName", "dGeneName", 
#' "dGeneAllele", "dFamilyTies", "dGeneNameTies", "dGeneAlleleTies", "jMaxResolved", 
#' "jFamilyName", "jGeneName", "jGeneAllele", "jFamilyTies", "jGeneNameTies", 
#' "jGeneAlleleTies", "vDeletion", "d5Deletion", "d3Deletion", "jDeletion", 
#' "n2Insertion", "n1Insertion", "vIndex", "n2Index", "dIndex", "n1Index", 
#' "jIndex", "estimatedNumberGenomes", "sequenceStatus", "cloneResolved", 
#' "vOrphon", "dOrphon", "jOrphon", "vFunction", "dFunction", "jFunction", 
#' "fractionNucleated".  
#'  
#' IMPORTANT: be aware that Adaptive has changed the 
#' column names of their files over time and if the headings of your files are 
#' inconsistent, then specify column = "all" or include all variations of the 
#' headings you want to important.  For example, column = c("count", 
#' "count (templates)", "count (reads)").  Also be aware that the "count" 
#' column previously reported the number of sequencing reads in earlier 
#' versions of ImmunoSEQ but now is equivalent to the 
#' "estimatedNumberGenomes" column.
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
library(tidyverse)

readAdaptiveV1 <- function(clone_file) {
    col_old = cols_only("nucleotide" = col_character(), "aminoAcid" = col_character(), 
        "count (reads)" = col_integer(), "frequencyCount (%)" = col_double(),
        "vGeneName" = col_character(), "dGeneName" = col_character(), "jGeneName" = col_character(), 
        "vFamilyName" = col_character(), "dFamilyName" = col_character(), "jFamilyName" = col_character(),
        "sequenceStatus" = col_character(), "estimatedNumberGenomes" = col_integer())
    col_std = c(nucleotide = "nucleotide", aminoAcid = "aminoAcid", count = "count (reads)", 
            frequencyCount = "frequencyCount (%)", vFamily = "vFamilyName", vGene = "vGeneName", 
            dFamily = "dFamilyName", dGene = "dGeneName", jFamily = "jFamilyName", jGene = "jGeneName",
            `function` = "sequenceStatus", estimatedNumberGenomes = "estimatedNumberGenomes")
    file_names <- tools::file_path_sans_ext(basename(clone_file))
    #names(file_list) <- file_names
    clone_frame <- readr::read_tsv(clone_file, col_types = col_old,
                                    na = c("", "NA", "Nan", "NaN"), trim_ws = TRUE,
                                    progress = show_progress())
    clone_frame <- clone_frame %>% rename(!!col_std)
    clone_frame <- clone_frame %>% mutate(`function`= str_replace(`function`, "In", "in-frame"))
    clone_frame <- clone_frame %>% add_column(sample=file_names)
    return(clone_frame)
}

readAdaptiveV2 <- function(clone_file) {
    col_old = cols_only("nucleotide" = col_character(), "aminoAcid" = col_character(), 
        "count (templates/reads)" = col_integer(), "frequencyCount (%)" = col_double(),
        "vGeneName" = col_character(), "dGeneName" = col_character(), "jGeneName" = col_character(), 
        "vFamilyName" = col_character(), "dFamilyName" = col_character(), "jFamilyName" = col_character(),
        "sequenceStatus" = col_character(), "estimatedNumberGenomes" = col_integer())
    col_std = c(nucleotide = "nucleotide", aminoAcid = "aminoAcid", count = "count (templates/reads)", 
            frequencyCount = "frequencyCount (%)", vFamily = "vFamilyName", vGene = "vGeneName", 
            dFamily = "dFamilyName", dGene = "dGeneName", jFamily = "jFamilyName", jGene = "jGeneName",
            `function` = "sequenceStatus", estimatedNumberGenomes = "estimatedNumberGenomes")
    file_names <- tools::file_path_sans_ext(basename(clone_file))
    #names(file_list) <- file_names
    clone_frame <- readr::read_tsv(clone_file, col_types = col_old,
                                    na = c("", "NA", "Nan", "NaN"), trim_ws = TRUE,
                                    progress = show_progress())
    clone_frame <- clone_frame %>% rename(!!col_std)
    clone_frame <- clone_frame %>% mutate(`function`= str_replace(`function`, "In", "in-frame"))
    clone_frame <- clone_frame %>% add_column(sample=file_names)
    return(clone_frame)
}

readAdaptiveV3 <- function(clone_file) {
    col_old = cols_only("nucleotide" = col_character(), "aminoAcid" = col_character(), 
        "count (templates)" = col_integer(), "frequencyCount (%)" = col_double(),
        "vGeneName" = col_character(), "dGeneName" = col_character(), "jGeneName" = col_character(), 
        "vFamilyName" = col_character(), "dFamilyName" = col_character(), "jFamilyName" = col_character(),
        "sequenceStatus" = col_character(), "estimatedNumberGenomes" = col_integer())
    col_std = c(nucleotide = "nucleotide", aminoAcid = "aminoAcid", count = "count (templates)", 
            frequencyCount = "frequencyCount (%)", vFamily = "vFamilyName", vGene = "vGeneName", 
            dFamily = "dFamilyName", dGene = "dGeneName", jFamily = "jFamilyName", jGene = "jGeneName",
            `function` = "sequenceStatus", estimatedNumberGenomes = "estimatedNumberGenomes")
    file_names <- tools::file_path_sans_ext(basename(clone_file))
    #names(file_list) <- file_names
    clone_frame <- readr::read_tsv(clone_file, col_types = col_old,
                                    na = c("", "NA", "Nan", "NaN"), trim_ws = TRUE,
                                    progress = show_progress())
    clone_frame <- clone_frame %>% rename(!!col_std)
    clone_frame <- clone_frame %>% mutate(`function`= str_replace(`function`, "In", "in-frame"))
    clone_frame <- clone_frame %>% add_column(sample=file_names)
    return(clone_frame)
}

readAdaptiveV4 <- function(clone_file) {
    col_old = cols_only("nucleotide.CDR3.in.lowercase." = col_character(), 
        "aminoAcid.CDR3.in.lowercase." = col_character(), "cloneCount" = col_integer(), 
        "clonefrequency...." = col_double(), "vGene" = col_character(), "dGene" = col_character(),
        "jGene" = col_character(), "fuction" = col_character())
    col_std = c(nucleotide = "nucleotide.CDR3.in.lowercase.", 
        aminoAcid = "aminoAcid.CDR3.in.lowercase.", count = "cloneCount", 
        frequencyCount = "clonefrequency....", vGene = "vGene", dGene = "dGene", jGene = "jGene",
        `function` = "fuction")
    file_names <- tools::file_path_sans_ext(basename(clone_file))
    clone_frame <- readr::read_tsv(clone_file, col_types = col_old,
                                    na = c("", "NA", "Nan", "NaN"), trim_ws = TRUE,
                                    progress = show_progress())
    clone_frame <- clone_frame %>% rename(!!col_std)
    clone_frame <- clone_frame %>% extract(nucleotide, c("nucleotide"), regex = "([acgt]+)")
    clone_frame <- clone_frame %>% mutate(nucleotide  = toupper(nucleotide))
    clone_frame <- clone_frame %>% extract(aminoAcid, c("aminoAcid"), regex = "([acdefghiklmnpqrstvwy]+)")
    clone_frame <- clone_frame %>% mutate(aminoAcid = toupper(aminoAcid))
    clone_frame <- clone_frame %>% separate(vGene, c("vFamily"), sep = '-', remove = FALSE)
    clone_frame <- clone_frame %>% separate(dGene, c("dFamily"), sep = '-', remove = FALSE)
    clone_frame <- clone_frame %>% separate(jGene, c("jFamily"), sep = '-', remove = FALSE)
    clone_frame <- clone_frame %>% add_column(sample=file_names)

    clone_frame <- clone_frame %>% group_by(nucleotide) %>% 
        summarize(aminoAcid = first(aminoAcid), `count` = sum(`count`), vGene = first(vGene), `function` = first(`function`),
        jGene = first(jGene), dGene = first(dGene)) %>% mutate(frequencyCount = `count` / sum(`count`))
    clone_frame <- clone_frame %>% add_colunn(estimatedNumberGenomes=`count`)
    return(clone_frame)
}

readBGIclone <- function(clone_file) {
    col_old <- cols_only("nucleotide(CDR3 in lowercase)" = col_character(), 
        "aminoAcid(CDR3 in lowercase)" = col_character(), "cloneCount" = col_integer(),
        "clonefrequency (%)" = col_double(), "vGene" = col_character(), "dGene" = col_character(),
        "jGene" = col_character(), "fuction" = col_character())
    col_std <- c(nucleotide = "nucleotide(CDR3 in lowercase)", 
        aminoAcid = "aminoAcid(CDR3 in lowercase)", count = "cloneCount", 
        frequencyCount = "clonefrequency (%)", vGene = "vGene", dGene = "dGene", jGene = "jGene",
        `function` = "fuction")
    file_names <- tools::file_path_sans_ext(basename(clone_file))
    clone_frame <- readr::read_tsv(clone_file, col_types = col_old,
                                    na = c("", "NA", "Nan", "NaN"), trim_ws = TRUE,
                                    progress = show_progress())
    clone_frame <- clone_frame %>% rename(!!col_std)
    clone_frame <- clone_frame %>% extract(nucleotide, c("nucleotide"), regex = "([acgt]+)")
    clone_frame <- clone_frame %>% mutate(nucleotide  = toupper(nucleotide))
    clone_frame <- clone_frame %>% extract(aminoAcid, c("aminoAcid"), regex = "([acdefghiklmnpqrstvwy]+)")
    clone_frame <- clone_frame %>% mutate(aminoAcid = toupper(aminoAcid))
    clone_frame <- clone_frame %>% separate(vGene, c("vFamily"), sep = '-', remove = FALSE)
    clone_frame <- clone_frame %>% separate(dGene, c("dFamily"), sep = '-', remove = FALSE)
    clone_frame <- clone_frame %>% separate(jGene, c("jFamily"), sep = '-', remove = FALSE)
    clone_frame <- clone_frame %>% add_column(sample=file_names)

    clone_frame <- clone_frame %>% group_by(nucleotide) %>% 
        summarize(aminoAcid = first(aminoAcid), `count` = sum(`count`), vGene = first(vGene), `function` = first(`function`),
        jGene = first(jGene), dGene = first(dGene)) %>% mutate(frequencyCount = `count` / sum(`count`))
    
    clone_frame <- clone_frame %>% add_colunn(estimatedNumberGenomes=`count`)
    return(clone_frame)
}

readImmunoSeq <- function(path, mode="adaptiveV2", recursive = FALSE) {
    file.paths1 <- list.files(path, full.names = TRUE, all.files = FALSE, 
                             recursive = recursive, pattern = ".tsv", 
                             include.dirs = FALSE)
    file.info <- file.info(file.paths1)
    file.paths2 <- rownames(file.info)[file.info$size > 0]
    if(!identical(file.paths1,file.paths2)){
        warning("One or more of the files you are trying to import has no sequences and will be ignored.", call. = FALSE)
    }

    if (mode == "adaptiveV2") {
        file.list <- suppressWarnings(plyr::llply(file.paths2, readAdaptiveV2, .progress = "text"))    
    } else if  (mode == "adaptiveV3") {
        file.list <- suppressWarnings(plyr::llply(file.paths2, readAdaptiveV3, .progress = "text"))    
    } else if (mode == "bgiClone") {
        file.list <- suppressWarnings(plyr::llply(file.paths2, readBGIclone, .progress = "text"))    
    } else if (mode == "adaptiveV1") {
        file.list <- suppressWarnings(plyr::llply(file.paths2, readAdaptiveV1, .progress = "text"))    
    }
    else if (mode == "adaptiveV4") {
        file.list <- suppressWarnings(plyr::llply(file.paths2, readAdaptiveV4, .progress = "text"))    
    }
    if(length(unique(plyr::llply(file.list, ncol))) > 1){
        warning("One or more of the files you are trying to import do not contain all the columns you specified.", call. = FALSE)
    }
    file.names <- gsub(".tsv", "", basename(file.paths2))
    names(file.list) <- file.names
    i <- 1
    for (i in 1:length(file.list)) {
        colnames(file.list[[i]]) <- ifelse(grepl("frequencyCount*", 
                                                 colnames(file.list[[i]])), 
                                           "frequencyCount", 
                                           colnames(file.list[[i]]))
        colnames(file.list[[i]]) <- ifelse(grepl("count*", 
                                                 colnames(file.list[[i]])), 
                                           "count", 
                                           colnames(file.list[[i]]))
    }
    return(file.list)
} 