% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/bhattacharyyaMatrix.R
\name{bhattacharyyaMatrix}
\alias{bhattacharyyaMatrix}
\title{Bhattacharyya matrix}
\usage{
bhattacharyyaMatrix(productive.seqs)
}
\arguments{
\item{productive.seqs}{A list data frames of productive sequences generated 
by the LymphoSeq function productiveSeq.  "frequencyCount" and "aminoAcid" 
are a required columns.}
}
\value{
A data frame of Bhattacharyya coefficients calculated from all 
pairwise comparisons from a list of sample data frames.  The Bhattacharyya 
coefficient is a measure of the amount of overlap between two samples.  The 
value ranges from 0 to 1 where 1 indicates the sequence frequencies are 
identical in the two samples and 0 indicates no shared frequencies.
}
\description{
Calculates the Bhattacharyya coefficient of all pairwise comparison from a 
list of data frames.
}
\examples{
file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")

file.list <- readImmunoSeq(path = file.path)

productive.aa <- productiveSeq(file.list, aggregate = "aminoAcid")

bhattacharyyaMatrix(productive.seqs = productive.aa)
}
\seealso{
\code{\link{pairwisePlot}} for plotting results as a heat map.
}
