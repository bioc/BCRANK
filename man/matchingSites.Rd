\name{matchingSites}
\docType{methods}
\alias{matchingSites}

\title{Report IUPAC consensus occurrences in a fasta file}

\description{
  This function reports all occurrences of a consensus sequence in a fasta
  file. It can be used to extract transcription factor binding sites
  predicted by BCRANK or other motif search methods.
}

\usage{
matchingSites(fafile, motifSequence, revComp=TRUE, strip.desc=TRUE)
}

\arguments{
  \item{fafile}{a ranked fasta file containing DNA sequences.} 
  \item{motifSequence}{a character vector in IUPAC coding representing a
	DNA sequence.}
  \item{revComp}{set to TRUE if the reverse complement also be matched.}
  \item{strip.desc}{An argument used by the function \code{\link[Biostrings:readFASTA]{readFASTA()}}.
	It specifies whether or not the ">" marking at the beginning of description lines should be removed.}
}

\value{
  Returns a data frame with positions, strand and DNA sequence for the
  matching sites.
}

\author{Adam Ameur, \email{adam.ameur@lcb.uu.se}}

\seealso{
  \code{\link{bcrank}}
}

\references{
 Ameur, A., Rada-Iglesias, A., Komorowski, J., Wadelius, C. Novel
 algorithm and ChIP-analysis identifies candidate functional SNPs. Submitted
}

\keyword{methods}

