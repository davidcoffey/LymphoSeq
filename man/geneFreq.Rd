% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/geneFreq.R
\name{geneFreq}
\alias{geneFreq}
\title{Gene frequencies}
\usage{
geneFreq(productive.nt, locus = "VDJ", family = FALSE)
}
\arguments{
\item{productive.nt}{A list of one or more data frames of productive sequences 
generated by the LymphoSeq function productiveSeq where the parameter 
aggregate is set to "nucleotide".}

\item{locus}{A character vector indicating which VDJ genes to include in the 
output.  Available options include "VDJ", "DJ", "VJ", "DJ", "V", "D", or "J".}

\item{family}{A Boolean value indicating whether or not family names instead 
of gene names are used.  If TRUE, then family names are used and if FALSE, 
gene names are used.}
}
\value{
Returns a data frame with the sample names, VDJ gene name, count, and 
\% frequency of the V, D, or J genes (each gene frequency should add to 
100\% for each sample).
}
\description{
Creates a data frame of VDJ gene counts and frequencies.
}
\examples{
file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")

file.list <- readImmunoSeq(path = file.path)

productive.nt <- productiveSeq(file.list = file.list, aggregate = "nucleotide")

geneFreq(productive.nt, locus = "VDJ", family = FALSE)

# Create a heat map of V gene usage
vfamilies <- geneFreq(productive.nt, locus = "V", family = TRUE)

require(reshape)

vfamilies <- reshape::cast(vfamilies, familyName ~ samples, value = "frequencyGene", sum)

rownames(vfamilies) <- as.character(vfamilies$familyName)

vfamilies$familyName <- NULL

RedBlue <- grDevices::colorRampPalette(rev(RColorBrewer::brewer.pal(11, "RdBu")))(256)

require(pheatmap)

pheatmap::pheatmap(vfamilies, color = RedBlue, scale = "row")

# Create a word cloud of V gene usage
vgenes <- geneFreq(productive.nt, locus = "V", family = FALSE)

require(wordcloud)

wordcloud::wordcloud(words = vgenes[vgenes$samples == "TRB_Unsorted_83", "geneName"], 
   freq = vgenes[vgenes$samples == "TRB_Unsorted_83", "frequencyGene"], 
	  colors = RedBlue)

# Create a cumulative frequency bar plot of V gene usage
vgenes <- geneFreq(productive.nt, locus = "V", family = FALSE)

require(ggplot2)

ggplot2::ggplot(vgenes, aes(x = samples, y = frequencyGene, fill = geneName)) +
  geom_bar(stat = "identity") +
  theme_minimal() + 
  scale_y_continuous(expand = c(0, 0)) + 
  guides(fill = guide_legend(ncol = 3)) +
  labs(y = "Frequency (\%)", x = "", fill = "") +
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust = 1))
}
