\name{bcrank}
\docType{methods}
\alias{bcrank}

\title{BCRANK: predicting binding site consensus from ranked DNA sequences}

\description{
  This function implements an algorithm for detection of short DNA
  sequences that are overrepresented in some part of the list. Starting
  from some initial consensus DNA sequence coded in IUPAC symbols, the
  method uses a heuristic search to improve the consensus until a local
  optimum is found. 
}

\usage{
bcrank(fafile, startguesses=c(), restarts=10, length=10,
       reorderings=500, silent=FALSE, plot.progress=FALSE,
       do.search=TRUE, use.P1=TRUE, use.P2=TRUE)
}

\arguments{
  \item{fafile}{a ranked fasta file containing DNA sequences.} 
  \item{startguesses}{a character vector with consensus sequences in
  IUPAC coding to be used as starting sequences in the search. If empty,
  random start guesses will be generated.}
  \item{restarts}{number restarts of the algorithm when using random start guesses.}
  \item{length}{legth of random start guess.}
  \item{reorderings}{number of random reorderings of the DNA sequences
	performed when calculating score.}
  \item{silent}{reports progress status if FALSE. }
  \item{plot.progress}{if TRUE, the progress is displayed in a plot.}
  \item{do.search}{if FALSE, no search is performed. In that case the
	start guesses are assigned with scores and reported as results.}
  \item{use.P1}{Use penalty for non-specific bases.}
  \item{use.P2}{Use penalty for motifs matching repetitive sequences}
}

\value{
  The method returns an objcet of class \code{\link[BCRANK:BCRANKresult-class]{BCRANKresult-class}}.
}

\author{Adam Ameur, \email{adam.ameur@lcb.uu.se}}

\seealso{
  \code{\link[BCRANK:BCRANKresult-class]{BCRANKresult-class}}
}

\examples{
## Load example fasta file  
fastaFile <- system.file("Exfiles/USF1_small.fa", package = "BCRANK") 
## Run BCRANK
\dontrun{BCRANKout <- bcrank(fastaFile, restarts=20)}
\dontshow{data(BCRANKout)}
## Show BCRANK results
toptable(BCRANKout)
## The top scoring result
topMotif <- toptable(BCRANKout,1)
## Plot BCRANK search path
plot(topMotif)
## Position Weight Matrix
pwm(topMotif, normalize=FALSE)
}

\references{
 Ameur, A., Rada-Iglesias, A., Komorowski, J., Wadelius, C. Novel
 algorithm and ChIP-analysis identifies candidate functional SNPs. Submitted
}

\keyword{methods}

