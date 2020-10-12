library(RColorBrewer)
library(data.table)
library(dplyr)
library(ggplot2)
library(gridExtra)
library(proxy)
library(ComplexHeatmap)
library(gridExtra)

# set the working directory and read the filenames in the TRA/KS directory into "ks.tra.filenames"
setwd("/Volumes/home/silo_warren/ImmunoSEQ/TRA/KS/")
getwd()
list.files(path = "/Volumes/home/silo_warren/ImmunoSEQ/TRA/KS/", pattern = ".tsv")
ks.tra.filenames <- list.files(path = "/Volumes/home/silo_warren/ImmunoSEQ/TRA/KS/", pattern = ".tsv")
ks.tra.files <- lapply(ks.tra.filenames, read.table, header = TRUE, sep = "\t", na.strings = c("", " ", "  ", "NA"))

# bind all of the data.frames in ks.tra.files into a single data.frame
all.ks.tra.seqs <- dplyr::bind_rows(ks.tra.files, .id = NULL)

# remove all of the rows with non-productive sequences (no CDR3 or containing a stop codon)
all.ks.tra.seqs.prod <- all.ks.tra.seqs[!is.na(all.ks.tra.seqs$aminoAcid),]
all.ks.tra.seqs.prod2 <- all.ks.tra.seqs.prod[!grepl("\\*", all.ks.tra.seqs.prod$aminoAcid),]

# define the data.frame output.tra.matrix and assign the unique TRA CDR3 sequences to the row names
# and TRA filenames to the columns (variables)
output.tra.matrix <- data.frame(row.names = unique(all.ks.tra.seqs.prod2$aminoAcid))

# initialize all positions in the output.tra.matrix
for (i in seq(1,length(ks.tra.filenames))){
  output.tra.matrix[,i] <- 0
}

#IF output.tra.matrix has already been generated saved - read it in HERE
output.tra.matrix <- data.table::fread("/Volumes/home/output.tra.matrix.csv", header = TRUE)

#assign the filenames to the variable (column) names in output.tra.matrix
colnames(output.tra.matrix) <- ks.tra.filenames

#shorten the column names in output.tra.matrix by taking only the first 7 characters
new.names <- substr(ks.tra.filenames, 1, 7)
colnames(output.tra.matrix) <- new.names

# define functions that will search for CDR3 sequences in ks.tra data.frames

findValInAminoAcid <- function(amino.string, list.of.dfs){
  #loop through and retrieve the full data frame if string is present, 0 otherwise
  locations <- list()
  for(i in seq(1,length(list.of.dfs),1)) {
    temp <- list.of.dfs[[i]]
    #subset the data frame to see if the value is present
    temp.sub<-temp[which(temp$aminoAcid==amino.string), ]
    if (nrow(temp.sub)==0){ #if the row length is 0 
      locations[[i]] <- 0
    }else{ # if it is present (row length > 0) add the full data frame
      locations[[i]] <- temp.sub
    }
  }
  return(locations)
}

findValInAminoAcid.binary <- function(locations){
  #loop through and return a binary list indicating the presence of the string
  locations.binary <- list()
  for (i in seq(1,length(locations),1)){
    temp <- locations[[i]]
    if (temp[1]==0){
      locations.binary[i]<-0
    }else{
      locations.binary[i]<-1
    }
  }
  return(locations.binary)
}

# take all of the unique TRA CDR3 sequences listed in output.tra.matrix and determine
#in which TRA datasets they are found
for(i in seq(1, nrow(output.tra.matrix), 1)){
  ans <- findValInAminoAcid(as.character(row.names(output.tra.matrix)[i]), ks.tra.files)
  loc.bin <- findValInAminoAcid.binary(ans)
  output.tra.matrix[i,] <- loc.bin
}

#read in the list of candidate KS public TRB sequences as KS.publicTRB.list
setwd("/User/ehwarren/")
getwd()
KS.publicTRB.list <- read.table(file = "/Users/ehwarren/KS.publicTRB.csv", header = TRUE, sep = ",",
                                colClasses = "character", na.strings = c("", " ", "  ", "NA"))
colnames(KS.publicTRB.list) <- c("aminoAcid", "mhcRestriction", "certainty", "specificity")

# define the trb filenames that will be loaded
ks.trb.filenames <- c("008_001_B.tsv", "008_001_D.tsv", "008_001_I.tsv", "008_002_A.tsv", "008_002_B.tsv", "008_002_C.tsv", "008_002_D.tsv",
                      "008_002_H.tsv", "008_002_I.tsv", "008_003_B.tsv", "008_003_3_D.tsv", "008_003_4_E.tsv", "008_003_5_F.tsv", "008_003_G.tsv",
                      "008_004_B.tsv", "008_004_C.tsv", "008_004_H.tsv", "008_004_I.tsv", "008_005_A.tsv", "008_005_B.tsv", "008_005_C.tsv",
                      "008_006_8_I.tsv", "008_006_C.tsv", "008_006_D.tsv", "008_006_H.tsv", "008_007_A.tsv", "008_007_H.tsv", "008_007_I.tsv",
                      "008_008_C.tsv", "008_008_D.tsv", "008_008_E.tsv", "008_008_F.tsv", "008_008_G.tsv", "008_008_H.tsv", "008_008_I.tsv",
                      "008_010_A_NAT.tsv", "008_010_B.tsv", "008_010_C.tsv", "008_011_9_J.tsv", "008_011_B.tsv", "008_011_D.tsv", "008_011_H.tsv",
                      "008_011_I.tsv", "008_012_A.tsv", "008_012_B.tsv", "008_012_C.tsv", "008_012_H.tsv", "008_012_I.tsv", "008_013_B.tsv",
                      "008_015_C.tsv", "008_015_E.tsv", "008_015_F.tsv", "008_015_H.tsv", "008_015_I.tsv", "008_016_A.tsv", "008_016_B.tsv",
                      "008_016_C.tsv", "008_016_D.tsv", "008_016_E.tsv", "008_016_F.tsv", "008_018_A_NAT.tsv", "008_018_B.tsv", "008_018_C.tsv",
                      "008_020_E.tsv", "008_021_A.tsv", "008_021_B.tsv", "008_021_C.tsv", "008_021_H.tsv", "008_021_I.tsv", "008_022_A.tsv",
                      "008_022_C.tsv", "008_023_A.tsv", "008_023_C.tsv", "008_023_D.tsv", "008_023_E.tsv", "008_024_B.tsv", "008_024_C.tsv",
                      "008_026_C.tsv", "008_028_C.tsv", "008_028_D.tsv", "008_028_I.tsv", "008_029_7_H.tsv", "008_029_B.tsv", "008_029_3_D.tsv",
                      "008_030_B.tsv", "008_030_3_D.tsv", "008_030_4_E.tsv", "008_030_5_F.tsv", "008_034_H.tsv", "008_034_I.tsv", "008_035_A.tsv",
                      "008_035_B.tsv", "008_035_C.tsv", "008_035_E.tsv", "008_035_F.tsv", "008_036_C.tsv", "008_037_7_H.tsv", "008_037_9_J.tsv",
                      "008_037_A_NAT.tsv", "008_037_B.tsv", "008_037_I.tsv", "008_039_A.tsv", "008_039_B.tsv", "008_039_C.tsv", "008_039_E.tsv",
                      "008_048_A.tsv", "008_048_C.tsv", "008_048_D.tsv", "008_048_H.tsv", "008_048_I.tsv", "008_049_A.tsv", "008_049_B.tsv",
                      "008_049_C.tsv", "008_052_B.tsv", "008_052_C.tsv", "008_052_F.tsv", "008_057_A.tsv", "008_057_B.tsv", "008_057_C.tsv",
                      "008_057_D.tsv", "008_059_A.tsv", "008_059_B.tsv", "008_059_D.tsv", "008_059_F.tsv", "008_060_C.tsv", "008_061_A.tsv",                           
                      "008_061_B.tsv", "008_061_C.tsv", "008_061_D.tsv", "008_061_E.tsv", "008_062_B.tsv", "008_062_C.tsv", "008_062_E.tsv",
                      "008_062_H.tsv", "008_062_I.tsv", "008_063_A.tsv", "008_063_E.tsv", "008_063_F.tsv", "008_066_H.tsv",  "008_066_I.tsv",
                      "008_067_B.tsv", "008_067_C.tsv", "008_067_D.tsv", "008_067_E.tsv", "008_067_F.tsv", "008_068_B_1.tsv", "008_068_C.tsv",
                      "008_069_B.tsv", "008_069_C.tsv", "008_001_A.tsv", "008_069_E.tsv", "008_069_F.tsv", "008_069_J.tsv", "008_069_NAT_A.tsv",
                      "008_070_A.tsv", "008_070_B.tsv", "008_070_D.tsv", "008_070_H.tsv", "008_071_B.tsv", "008_071_C.tsv", "008_071_H.tsv",
                      "008_071_I.tsv", "008_074_B.tsv", "008_074_C.tsv", "008_074_D.tsv", "008_080_A.tsv", "008_080_B.tsv", "008_080_C.tsv",
                      "008_085_B.tsv", "008_088_D.tsv", "008_091_B.tsv", "008_092_D.tsv", "008_092_H.tsv", "008_093_B.tsv", "008_093_C.tsv",
                      "008_095_A.tsv", "008_095_B.tsv", "008_098_A.tsv", "008_098_B.tsv", "008_098_C.tsv", "008_099_D.tsv", "008_101_D.tsv",
                      "008_115_B.tsv", "008_003_A_NAT.tsv", "008_115_D.tsv", "008_122_B.tsv", "008_124_A.tsv", "008_124_B.tsv", "008_124_D.tsv",
                      "008_127_A.tsv", "008_127_C.tsv", "008_127_G.tsv", "008_127_H.tsv", "008_130_A.tsv", "008_130_B.tsv", "008_130_E.tsv",
                      "008_130_F.tsv", "008_138_B.tsv", "008_138_C.tsv", "008_138_F.tsv", "008_138_G.tsv", "008_140_B.tsv", "008_141_H.tsv",
                      "008_143_A.tsv", "008_143_B.tsv", "008_143_C.tsv", "008_143_D.tsv", "008_145_A.tsv", "008_145_B.tsv", "008_145_C.tsv",
                      "008_146_A.tsv", "008_146_B.tsv", "008_146_C.tsv", "008_146_H.tsv", "008_146_I.tsv", "008_149_C.tsv", "008_149_F.tsv", 
                      "008_149_G.tsv", "008_156_A.tsv", "008_156_B.tsv", "008_156_C.tsv", "008_156_G.tsv", "008_156_H.tsv", "008_166_A.tsv",
                      "008_166_B.tsv", "008_166_C.tsv", "008_166_G.tsv", "008_166_H.tsv", "008_175_B.tsv", "008_175_C.tsv", "008_175_D.tsv",
                      "008_175_E.tsv", "008_175_F.tsv", "008_175_G.tsv", "008_176_A.tsv", "008_176_B.tsv", "008_176_C.tsv", "008_176_D.tsv",                          
                      "008_178_A.tsv", "008_178_C.tsv", "008_180_C.tsv", "008_187_B.tsv", "008_187_C.tsv", "008_187_F.tsv", "008_187_G.tsv",
                      "008_188_A.tsv", "008_189_A.tsv", "008_189_B.tsv", "008_192_A.tsv", "008_192_D.tsv", "008_196_B.tsv", "008_196_C.tsv",                          
                      "008_196_D.tsv", "008_196_F.tsv", "008_197_D.tsv", "008_208_A_NAT.tsv", "008_208_D_TRB.tsv", "008_208_KS_SingleCell_BulkSuspension.tsv",
                      "008_210_A_NAT.tsv", "008_210_E_TRB.tsv", "008_211_A_TRB.tsv", "008_211_D_TRB.tsv", "008_213_A_TRB.tsv", "008_213_D_TRB.tsv",                          
                      "008_213_KS_SingleCell_Suspension.tsv", "008_215_A_TRB.tsv", "008_215_D_TRB.tsv", "008_216_A_TRB.tsv", "008_216_D_TRB.tsv", "008_217_A_TRB.tsv",                       
                      "008_217_D_TRB.tsv", "008_217_PBMC_TRB.tsv", "008_220_A_TRB.tsv", "008_220_D_TRB.tsv", "008_221_A_TRB.tsv", "008_221_D_TRB.tsv", "008_225_A_TRB.tsv")

#define the data.frame output.trb.matrix and assign the candidate public TRB CDR3 sequences to the row names
output.trb.matrix <- data.frame(row.names = KS.publicTRB.list$aminoAcid, check.rows = TRUE, check.names = TRUE)

#initialize all positions in the output.trb.matrix
for (i in seq(1,length(ks.trb.filenames))){
  output.trb.matrix[,i] <- 0
}

#assign the filenames to the variable (column) names in output.trb.matrix
colnames(output.trb.matrix) <- ks.trb.filenames

#shorten the column names in output.trb.matrix by taking only the first 7 characters
new.names <- substr(ks.trb.filenames, 1, 7)
colnames(output.trb.matrix) <- new.names

#set the working directory and read the files in "ks.trb.filenames" into ks.trb.files
setwd("/Volumes/home/silo_warren/ImmunoSEQ/TRB/KS/")
getwd()
ks.trb.files <- lapply(ks.trb.filenames, read.table, header = TRUE, sep = "\t", na.strings = c("", " ", "  ", "NA"))

#determine which TRB datasets the candidate public TRB CDR3 sequences are found in
for(i in seq(1, nrow(output.trb.matrix), 1)){
  ans <- findValInAminoAcid(as.character(row.names(output.trb.matrix)[i]), ks.trb.files)
  loc.bin <- findValInAminoAcid.binary(ans)
  output.trb.matrix[i,] <- loc.bin
}

pheatmap(output.trb.matrix, color = cm.colors(256), scale = "none", cluster_rows = FALSE, cluster_cols = FALSE,
         legend = FALSE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 14, fontsize_row = 6, fontsize_col = 6, border_color = NA,
         main = "Selected CDR3 Sequences in All KS TRB Samples, No Clustering")

pheatmap(output.trb.matrix, color = cm.colors(256), scale = "none", cluster_rows = TRUE, cluster_cols = FALSE,
         legend = FALSE, show_rownames = TRUE, show_colnames = TRUE, fontsize = 14, fontsize_row = 6, fontsize_col = 6, border_color = NA,
         main = "Selected CDR3 Sequences in All KS TRB Samples, Clustering by Sequence")

# take a specific row of output.trb.matrix and save it as test.vector
test.vector <- as.matrix(output.trb.matrix[65,])

# compute the Jaccard distance between test.vector and every row of output.tra.matrix
similar.results <- dist(as.matrix(output.tra.matrix), y = test.vector, method = "Jaccard", diag = FALSE, upper = FALSE,
                        pairwise = FALSE, by_rows = TRUE, auto_convert_data_frames = TRUE)
similar.results <- as.matrix(similar.results)

# find the row with the closest Jaccard distance to the test vector
min(similar.results)
sel = which(min(similar.results[,1]) == similar.results[,1])
sel
output.tra.matrix[sel,]


for(i in seq(1, nrow(output.trb.matrix), 1)){
  test.vector <- as.matrix(output.trb.matrix[i,])
  similar.results <- dist(as.matrix(output.tra.matrix), y = test.vector, method = "Jaccard", diag = FALSE, upper = FALSE,
                          pairwise = FALSE, by_rows = TRUE, auto_convert_data_frames = TRUE)
  sel = which(min(similar.results[,1]) == similar.results[,1])
  print(sel)
}

ht_opt("heatmap_column_names_gp" = gpar(fontsize = 5), "heatmap_row_names_gp" = gpar(fontsize = 20, fontface = "bold"))

ht1 <- Heatmap(as.matrix(output.trb.matrix[37,]), col = cm.colors(256), cluster_rows = FALSE, cluster_columns = FALSE,
               column_title = "Finding best TRA match for TRB sequence CASSIAGHEQYF", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               rect_gp = gpar(col = "black", lty = 1, lwd = 1))
ht2 <- Heatmap(as.matrix(output.tra.matrix[7737,]), col = cm.colors(256), cluster_rows = FALSE, cluster_columns = FALSE,
               rect_gp = gpar(col = "black", lty = 1, lwd = 1))

ht_list <- ht1 %v% ht2
draw(ht_list, show_heatmap_legend = FALSE, gap = unit(5, "mm"))

ht3 <- Heatmap(as.matrix(output.trb.matrix[65,]), col = cm.colors(256), cluster_rows = FALSE, cluster_columns = FALSE,
               column_title = "Finding best TRA match for TRB sequence CASSLWGGPSNEQFF", column_title_gp = gpar(fontsize = 20, fontface = "bold"),
               rect_gp = gpar(col = "black", lty = 1, lwd = 1))
ht4 <- Heatmap(as.matrix(output.tra.matrix[4249,]), col = cm.colors(256), cluster_rows = FALSE, cluster_columns = FALSE,
               rect_gp = gpar(col = "black", lty = 1, lwd = 1))

ht_list <- ht3 %v% ht4
draw(ht_list, show_heatmap_legend = FALSE, gap = unit(5, "mm"))


best.match <- output.tra.matrix[grepl("CAYSGAGSYQLTF", row.names(output.tra.matrix)),]


