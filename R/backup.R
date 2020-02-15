library(tidyverse)
library(proxy)
library(Matrix)

# Define custom functions
tic()
# Function to identify the different adaptive output file formats
getFileType <- function(clone_file) {
    adaptiveV1 <- c("nucleotide", "aminoAcid", "count (reads)", "frequencyCount (%)", "vGeneName", 
                    "dGeneName", "jGeneName")
    
    adaptiveV2 <- c("nucleotide", "aminoAcid", "count (templates/reads)", "frequencyCount (%)", "vGeneName", 
                    "dGeneName", "jGeneName")
    
    adaptiveV3 <- c("nucleotide", "aminoAcid", "count (templates)", "frequencyCount (%)", "vGeneName", 
                    "dGeneName", "jGeneName")
    
    columns <- invisible(colnames(readr::read_tsv(clone_file, n_max=1, col_types = cols())))
    if (all(adaptiveV1 %in% columns)) {
        file_type <- "adaptiveV1"
        header_list <- cols_only(nucleotide = "c", aminoAcid = "c", `count (reads)` = "i", 
                                 `frequencyCount (%)` = "d", vGeneName = "c", dGeneName = "c", 
                                 jGeneName = "c")
    } else if (all(adaptiveV2 %in% columns)) {
        file_type <- "adaptiveV2"
        header_list <- cols_only(nucleotide = "c", aminoAcid = "c", `count (templates/reads)` = "i", 
                                 `frequencyCount (%)` = "d", vGeneName = "c", dGeneName = "c", 
                                 jGeneName = "c")
    } else if (all(adaptiveV3 %in% columns)) {
        file_type <- "adaptiveV3"
        header_list <- cols_only(nucleotide = "c", aminoAcid = "c", `count (templates)` = "i", 
                                 `frequencyCount (%)` = "d", vGeneName = "c", dGeneName = "c", 
                                 jGeneName = "c")
    }
    ret_val <- list(file_type, header_list)
    return(ret_val)
}

# Fucntion to rename columns to standardize the nomeclature 
getStandard <- function(file_type) {
    adaptiveV1 <- c(count = "count (reads)", frequencyCount = "frequencyCount (%)")
    
    adaptiveV2 <- c(count = "count (templates/reads)", frequencyCount = "frequencyCount (%)")
    
    adaptiveV3 <- c(count = "count (templates)", frequencyCount = "frequencyCount (%)")
    

    type_hash <- list("adaptiveV1"=adaptiveV1, "adaptiveV2"=adaptiveV2, "adaptiveV3"=adaptiveV3)
    return(type_hash[[file_type]])
}

# Read adaptive files into tibbles, select only required columns, rename columns and add sample name column
readFiles <- function(clone_file) {
    file_info <- getFileType(clone_file)
    file_type <- file_info[[1]]
    header_list <- file_info[[2]]
    col_std <- getStandard(file_type)
    col_old = header_list
    file_names <- tools::file_path_sans_ext(basename(clone_file)) %>% str_split("_TRA|_TRB", n = 2, simplify = TRUE)
    clone_frame <- readr::read_tsv(clone_file, col_types = col_old,
                                   na = c("", "NA", "Nan", "NaN"), trim_ws = TRUE)
    clone_frame <- clone_frame %>% dplyr::rename(!!!col_std)
    clone_frame <- clone_frame %>% group_by(aminoAcid) %>%  summarize( `count` = sum(`count`), vGeneName = first(vGeneName),
                      jGeneName = first(jGeneName), dGeneName = first(dGeneName)) %>% mutate(frequencyCount = `count` / sum(`count`))
    clone_frame <- clone_frame %>% add_column(sample=file_names[,1])
    return(clone_frame)
}

# Get a list of ImmunoSeq files from TRA and TRB sequencing
getTCRMatrix <- function(tra_path, trb_path) {

    ks.tra.filenames <- list.files(path = tra_path, pattern = ".tsv", full.names = TRUE)
    ks.trb.filenames <- list.files(path = trb_path, pattern = ".tsv", full.names = TRUE)
    # Sample file names are formatted such that files from TRA and TRB folder have the same name
    ks.tra.samplename <- tools::file_path_sans_ext(basename(ks.tra.filenames)) %>% str_split("_TRA", n = 2, simplify = TRUE) %>% as_tibble() %>% select(V1)
    ks.trb.samplename <- tools::file_path_sans_ext(basename(ks.trb.filenames)) %>% str_split("_TRB", n = 2, simplify = TRUE) %>% as_tibble() %>% select(V1)
    # Create a named vector with the keys being the formatted sample names
    names(ks.tra.filenames) <- ks.tra.samplename$V1
    names(ks.trb.filenames) <- ks.trb.samplename$V1
    # Get a list of sample names for which TRA and TRB sequencing files are available
    ks.tra.trb.matched <- inner_join(ks.tra.samplename, ks.trb.samplename, by="V1")
    # Create list of files for which both TRA and TRB sequencing files are avaiable
    ks.tra.filenames.matched <- ks.tra.filenames[ks.tra.trb.matched$V1]
    ks.trb.filenames.matched <- ks.trb.filenames[ks.tra.trb.matched$V1]
    # Read all matched TRA and TRB files into two separate tibbles
    all.ks.tra.files <- ks.tra.filenames.matched %>% map(readFiles) %>% reduce(rbind)
    all.ks.trb.files <- ks.trb.filenames.matched %>% map(readFiles) %>% reduce(rbind)
    # Select productive sequences and filter entries on length of aminoAcid seqeunces, and count. Add a column called productive whose value is set to 1 for all productive sequences
    all.ks.tra.prod <- all.ks.tra.files %>% filter(!grepl("\\*", aminoAcid) & str_length(aminoAcid) > 7) %>% add_column(productive = 1) 
    all.ks.trb.prod <- all.ks.trb.files %>% filter(!grepl("\\*", aminoAcid) & str_length(aminoAcid) > 7)%>% add_column(productive = 1)
    # Pivot the tibble on aminoAcid and sample columns from the two tibbles, productive column is used to fill values
    # The resultant matrix is a binary matrix with 1 if a productive aminoAcid sequences is in a sample anc 0 otherwise
    ks.tra.matrix <- all.ks.tra.prod %>% pivot_wider(id_cols = aminoAcid, names_from = sample, values_from = productive, values_fill = list(productive =0))
    ks.trb.matrix <- all.ks.trb.prod %>% pivot_wider(id_cols = aminoAcid, names_from = sample, values_from = productive, values_fill = list(productive =0)) 
    # Filtering the tibbles by count and aminoAcid length results in some samples being omitted from the analysis
    # In order to calculate Jaccard distance between TRA and TRB sequenes we need to rematch the samples from the TRA and TRB matrices
    ks.tra.colnames <- colnames(ks.tra.matrix) %>% as_tibble()
    ks.trb.colnames <- colnames(ks.trb.matrix) %>% as_tibble()
    
    ks.tra.trb.matched <- inner_join(ks.tra.colnames, ks.trb.colnames, by="value")
    
    ks.tra.matrix <- ks.tra.matrix %>% select(ks.tra.trb.matched$value)
    ks.trb.matrix <- ks.trb.matrix %>% select(ks.tra.trb.matched$value)
    # Transform tibble to sparse matrix and tranpose TRB matrix for matrix multiplication
    tra.matrix <- ks.tra.matrix %>% select(-aminoAcid) %>% as.matrix()
    #tra.matrix <- Matrix(tra.matrix, sparse = T)
    trb.matrix <- ks.trb.matrix %>% select(-aminoAcid) %>% as.matrix() 
    #trb.matrix <- t(Matrix(trb.matrix, sparse = T))
    tra.list <- ks.tra.matrix$aminoAcid
    trb.list <- ks.trb.matrix$aminoAcid
    tcr.list <- list(tra.matrix, trb.matrix, tra.list, trb.list)
    return(tcr.list)
}

ks_tcr_list <- getTCRMatrix("/Volumes/home/fuser/kstme/immunoSeq/rawData/tra", "/Volumes/home/fuser/kstme/immunoSeq/rawData/trb")
ks_tra_mtx <- ks_tcr_list[[1]]
ks_trb_mtx <- ks_tcr_list[[2]]
ks_tra_list <- ks_tcr_list[[3]]
ks_trb_list <- ks_tcr_list[[4]]

ks_paired <- ks_tra %*% t(ks_trb)
colnames(ks_paired) <- ks_trb_list
ks_paired <- as_tibble(ks_paired)
ks_paired <- ks_paired %>% add_column(tra_sequence = ks_tra_list) %>% select(tra_sequence, everything())
ks_paired <- ks_paired %>% mutate_all(na_if, 1)
ks_paired <- ks_paired %>% pivot_longer(-tra_sequence, names_to = "trb_sequence", values_to = "jaccard", values_drop_na = TRUE)
ks_paired_filtered <- ks_paired %>% filter(jaccard > 2)

readClonotypes <- function(clone_file){
    clone_frame <- read_csv(clone_file)
    sample = tools::file_path_sans_ext(basename(clone_file))
    clone_frame <- clone_frame %>% add_column(sample = sample)
    return(clone_frame)
}

clonotypes <- list.files("/Volumes/home/fuser/kstme/clonotypes", pattern = ".csv", full.names = TRUE) %>% map(readClonotypes) %>% reduce(rbind)
clonotypes <- clonotypes %>% separate_rows(cdr3s_aa, cdr3s_nt, sep = ";")
clonotypes <- clonotypes %>% separate(cdr3s_aa, into = c("sequence_type", "aminoAcid"), sep = ":") %>% separate(cdr3s_nt, into = c("nuc_type", "nucleotide")) %>% select(-nuc_type)
clonotypes_tra <- clonotypes %>% filter(sequence_type == 'TRA')
clonotypes_trb <- clonotypes %>% filter(sequence_type == 'TRB')
clonotypes_paired <- full_join(clonotypes_tra, clonotypes_trb, by="clonotype_id")
clonotypes_paired <- clonotypes_paired %>% select(aminoAcid.x, aminoAcid.y, frequency.x, frequency.y) %>% rename(tra_sequence = "aminoAcid.x", trb_sequence = "aminoAcid.y", tra_count = "frequency.x", trb_count = "frequency.y")

adaptive_v_10x <- left_join(clonotypes_paired, ks_paired, by=c("tra_sequence", "trb_sequence"))

