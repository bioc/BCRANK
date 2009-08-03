\name{BCRANKresult-class}
\docType{class}
\alias{BCRANKresult-class}
\alias{BCRANKresult}
\alias{toptable}
\alias{fname}

\title{Class "BCRANKresult"}

\description{Holds the results from running
  \code{\link[BCRANK:bcrank]{bcrank}. Contains a number of
	\code{\link[BCRANK:BCRANKsearch-class]{BCRANKsearch-class}} object,
	one for each restart of the bcrank search.}
}

\section{Slots}{
  \describe{
    \item{\code{fname}:}{the name of the fasta file used for running bcrank.}
    \item{\code{toplist}:}{a list of \code{BCRANKsearch-class} objects,
	  ranked by their scores.}
	\item{\code{funCall}:}{the function call that was made to bcrank.}
	\item{\code{nrSeqs}:}{number of sequences in the fasta input file.}
	\item{\code{restarts}:}{number of restarts used in the bcrank
	  search.}
  }
}

\section{Methods}{
  \describe{
    \item{fname}{\code{signature(object = "BCRANKmatch")}: Returns the
	  fasta file name.}
	\item{toplist}{\code{signature(object = "BCRANKmatch", i=NULL)}: If
	  \code{i} is NULL, returns a data frame containing consensus and
	  score for the results for each restart of the bcrank
	  search. Otherwise, the i'th \code{BCRANKsearch-class} object in
	  the toplist is returned.}
  }
}

\author{Adam Ameur, \email{adam.ameur@genpat.uu.se}}

\seealso{
  \code{\link[BCRANK:bcrank]{bcrank}},
  \code{\link[BCRANK:BCRANKsearch-class]{BCRANKsearch-class}},
}

\keyword{classes}
