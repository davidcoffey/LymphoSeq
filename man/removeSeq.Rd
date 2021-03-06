% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/removeSeq.R
\name{removeSeq}
\alias{removeSeq}
\title{Remove sequence}
\usage{
removeSeq(file.list, sequence)
}
\arguments{
\item{file.list}{A list of data frames imported using the LymphoSeq function 
readImmunoSeq.  "aminoAcid", "count", and "frequencyCount" are required columns.}

\item{sequence}{A character vector of one or more amino acid sequences to 
remove from the list of data frames.}
}
\value{
Returns a list of data frames like the one imported except all rows 
with the specified amino acid sequence are removed.  The frequencyCount is 
recalculated.
}
\description{
Removes an amino acid sequence and associated data from all instances within 
a list of data frames and then recomputes the frequencyCount.
}
\examples{
file.path <- system.file("extdata", "TCRB_sequencing", package = "LymphoSeq")

file.list <- readImmunoSeq(path = file.path)

searchSeq(list = file.list, sequence = "CASSDLIGNGKLFF")

cleansed <- removeSeq(file.list = file.list, sequence = "CASSDLIGNGKLFF")

searchSeq(list = cleansed, sequence = "CASSDLIGNGKLFF")
}
