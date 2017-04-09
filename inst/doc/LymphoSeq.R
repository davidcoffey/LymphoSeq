## ---- results = "hide"---------------------------------------------------
library(LymphoSeq)

TCRB.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")

TCRB.list <- readImmunoSeq(path = TCRB.path)

## ---- comment = ""-------------------------------------------------------
names(TCRB.list)

## ---- comment = ""-------------------------------------------------------
lapply(TCRB.list, dim)

## ---- comment = ""-------------------------------------------------------
CMV <- TCRB.list[grep("CMV", names(TCRB.list))]
names(CMV)
TRB_Unsorted_0 <- TCRB.list[["TRB_Unsorted_0"]]
head(TRB_Unsorted_0)

## ---- comment = ""-------------------------------------------------------
TCRB.metadata <- read.csv(system.file("extdata", "TCRB_metadata.csv", package = "LymphoSeq"))
TCRB.metadata
selected <- as.character(TCRB.metadata[TCRB.metadata$phenotype == "Unsorted" & 
                                 TCRB.metadata$day > 300, "samples"])
TCRB.list.selected <- TCRB.list[selected]
names(TCRB.list.selected)

## ---- results = "hide"---------------------------------------------------
productive.TRB.aa <- productiveSeq(file.list = TCRB.list, aggregate = "aminoAcid", 
                               prevalence = FALSE)

## ---- results = "hide"---------------------------------------------------
productive.TRB.nt <- productiveSeq(file.list = TCRB.list, aggregate = "nucleotide", 
                               prevalence = FALSE)

## ---- comment = ""-------------------------------------------------------
head(TCRB.list[["TRB_Unsorted_0"]])

## ---- comment = ""-------------------------------------------------------
head(productive.TRB.aa[["TRB_Unsorted_0"]])

## ---- comment = ""-------------------------------------------------------
head(productive.TRB.nt[["TRB_Unsorted_0"]])

## ---- comment = ""-------------------------------------------------------
clonality(file.list = TCRB.list)

## ---- results = "hide"---------------------------------------------------
IGH.path <- system.file("extdata", "IGH_sequencing", package = "LymphoSeq")

IGH.list <- readImmunoSeq(path = IGH.path)

## ---- comment = ""-------------------------------------------------------
clonalRelatedness(list = IGH.list, editDistance = 10)

## ---- results = "hide"---------------------------------------------------
productive.IGH.nt <- productiveSeq(file.list = IGH.list, aggregate = "nucleotide")

## ---- fig.width = 7, fig.height = 8, comment = ""------------------------
phyloTree(list = productive.IGH.nt, sample = "IGH_MVQ92552A_BL", type = "nucleotide", 
         layout = "rectangular")

## ---- comment = ""-------------------------------------------------------
alignSeq(list = productive.IGH.nt, sample = "IGH_MVQ92552A_BL", type = "aminoAcid", 
         method = "ClustalW", output = "consule")

## ---- comment = ""-------------------------------------------------------
searchSeq(list = productive.TRB.aa, sequence = "CASSPVSNEQFF", type = "aminoAcid", 
          match = "global", editDistance = 0)

## ---- comment = ""-------------------------------------------------------
published <- searchPublished(list = productive.TRB.aa)
head(published)

## ---- fig.width = 7, fig.height = 7, comment = ""------------------------
lorenzCurve(samples = names(productive.TRB.aa), list = productive.TRB.aa)

## ---- fig.width = 7, fig.height = 5, comment = ""------------------------
topSeqsPlot(list = productive.TRB.aa, top = 10)

## ---- comment = ""-------------------------------------------------------
bhattacharyya.matrix <- bhattacharyyaMatrix(productive.seqs = productive.TRB.aa)
bhattacharyya.matrix[,1:2]
similarity.matrix <- similarityMatrix(productive.seqs = productive.TRB.aa)
similarity.matrix[,1:2]

## ---- fig.width = 6.5, fig.height = 5, comment = ""----------------------
pairwisePlot(matrix = bhattacharyya.matrix)

## ---- comment = ""-------------------------------------------------------
common <- commonSeqs(samples = c("TRB_Unsorted_0", "TRB_Unsorted_32"), 
                    productive.aa = productive.TRB.aa)
head(common)

## ---- fig.width = 4, fig.height = 4, comment = ""------------------------
commonSeqsVenn(samples = c("TRB_Unsorted_32", "TRB_Unsorted_83"), 
               productive.seqs = productive.TRB.aa)

## ---- fig.width = 4, fig.height = 4, comment = ""------------------------
commonSeqsVenn(samples = c("TRB_Unsorted_0", "TRB_Unsorted_32", "TRB_Unsorted_83"), 
               productive.seqs = productive.TRB.aa)

## ---- fig.width = 4, fig.height = 4, comment = ""------------------------
commonSeqsPlot("TRB_Unsorted_32", "TRB_Unsorted_83", 
               productive.aa = productive.TRB.aa)

## ---- fig.width = 7, fig.height = 5, comment = ""------------------------
commonSeqsBar(productive.aa = productive.TRB.aa, 
              samples = c("TRB_CD4_949", "TRB_CD8_949", 
                          "TRB_Unsorted_949", "TRB_Unsorted_1320"), 
              color.sample = "TRB_CD8_949")

## ---- comment = "", warning = FALSE--------------------------------------
differentialAbundance(list = productive.TRB.aa, 
                      sample1 = "TRB_Unsorted_949", 
                      sample2 = "TRB_Unsorted_1320", 
                      type = "aminoAcid", q = 0.01)

## ---- comment = ""-------------------------------------------------------
unique.seqs <- uniqueSeqs(productive.aa = productive.TRB.aa)
head(unique.seqs)
sequence.matrix <- seqMatrix(productive.aa = productive.TRB.aa, sequences = unique.seqs$aminoAcid)
head(sequence.matrix)[1:6]

## ---- comment = ""-------------------------------------------------------
top.freq <- topFreq(productive.aa = productive.TRB.aa, percent = 0.1)
head(top.freq)

## ---- comment = ""-------------------------------------------------------
top.freq <- topFreq(productive.aa = productive.TRB.aa, percent = 0)
top.freq.matrix <- merge(top.freq, sequence.matrix)
head(top.freq.matrix)[1:12]

## ---- fig.width = 7, fig.height = 5--------------------------------------
cloneTrack(sequence.matrix = sequence.matrix, 
           productive.aa = productive.TRB.aa, 
           map = c("TRB_CD4_949", "TRB_CD8_949"), 
           label = c("CD4", "CD8"), 
           track = "CASSPPTGERDTQYF", 
           unassigned = FALSE)

## ---- comment = ""-------------------------------------------------------
vGenes <- geneFreq(productive.nt = productive.TRB.nt, locus = "V", family = TRUE)
head(vGenes)

## ---- fig.width = 4, fig.height = 4--------------------------------------
top.seqs <- topSeqs(productive.seqs = productive.TRB.nt, top = 1)
chordDiagramVDJ(sample = top.seqs, 
                association = "VJ", 
                colors = c("darkred", "navyblue"))

## ---- fig.width = 4, fig.height = 4, warning = FALSE, message = FALSE----
vGenes <- geneFreq(productive.nt = productive.TRB.nt, locus = "V", family = TRUE)
library(RColorBrewer)
library(grDevices)
RedBlue <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)
library(wordcloud)
wordcloud::wordcloud(words = vGenes[vGenes$samples == "TRB_Unsorted_83", "familyName"], 
                     freq = vGenes[vGenes$samples == "TRB_Unsorted_83", "frequencyGene"], 
                     colors = RedBlue)

## ---- fig.width = 5, fig.height = 7, warning = FALSE, message = FALSE----
library(reshape)
vGenes <- reshape::cast(vGenes, familyName ~ samples, value = "frequencyGene", sum)
rownames(vGenes) = as.character(vGenes$familyName)
vGenes$familyName = NULL
library(pheatmap)
pheatmap::pheatmap(vGenes, color = RedBlue, scale = "row")

## ---- fig.width = 7, fig.height = 5.6, warning = FALSE, message = FALSE----
vGenes <- geneFreq(productive.nt = productive.TRB.nt, locus = "V", family = TRUE)
library(ggplot2)
multicolors <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(9, "Set1")))(28)
ggplot2::ggplot(vGenes, aes(x = samples, y = frequencyGene, fill = familyName)) +
  geom_bar(stat = "identity") +
  theme_minimal() + 
  scale_y_continuous(expand = c(0, 0)) + 
  guides(fill = guide_legend(ncol = 2)) +
  scale_fill_manual(values = multicolors) + 
  labs(y = "Frequency (%)", x = "", fill = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))

## ---- comment = ""-------------------------------------------------------
searchSeq(list = productive.TRB.aa, sequence = "CASSDLIGNGKLFF")
cleansed <- removeSeq(file.list = productive.TRB.aa, sequence = "CASSDLIGNGKLFF")
searchSeq(list = cleansed, sequence = "CASSDLIGNGKLFF")

## ------------------------------------------------------------------------
TRB_949_Merged <- mergeFiles(samples = c("TRB_CD4_949", "TRB_CD8_949"), 
                                file.list = TCRB.list)

## ------------------------------------------------------------------------
sessionInfo()

