\name{BCRANKmatch-class}
\docType{class}
\alias{BCRANKmatch-class}
\alias{BCRANKmatch}
\alias{consensus}
\alias{bcrankScore}
\alias{matchVector}

\title{Class "BCRANKmatch"}

\description{Holds the \code{\link[BCRANK:bcrank]{bcrank}} score for one
	IUPAC consensus sequence. Several objects of this class are
	collected in a
	\code{\link[BCRANK:BCRANKsearch-class]{BCRANKsearch-class}} object}
\section{Objects from the Class}{
  Objects are not intended to be created directly but as a result from
  running \code{\link[BCRANK:bcrank]{bcrank}}.
}

\section{Slots}{
  \describe{
    \item{\code{consensus}:}{consensus sequence in IUPAC coding}
    \item{\code{bcrankScore}:}{bcrank  score for the consensus}
  }
}

\section{Methods}{
  \describe{
    \item{consensus}{\code{signature(object = "BCRANKmatch")}: Returns
	  the consensus sequence.}
	\item{bcrankScore}{\code{signature(object = "BCRANKmatch")}: Returns
	  the bcrank score.}
	\item{matchVector}{\code{signature(object = "BCRANKmatch")}: Returns
	  a vector with 0's (no match) and 1's (match) of same length and
	  order as the ranked DNA sequences.}
  }
}

\author{Adam Ameur, \email{adam.ameur@lcb.uu.se}}

\seealso{
  \code{\link[BCRANK:bcrank]{bcrank}},
  \code{\link[BCRANK:BCRANKsearch-class]{BCRANKsearch-class}}
}

\keyword{classes}
