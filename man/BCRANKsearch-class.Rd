\name{BCRANKsearch-class}
\docType{class}
\alias{BCRANKsearch-class}
\alias{BCRANKsearch}
\alias{searchPath}
\alias{pwm}
\alias{plot,BCRANKsearch,missing-method}

\title{Class "BCRANKsearch"}

\description{Holds the whole search path from a single
  \code{\link[BCRANK:bcrank]{bcrank}} run. Each individual search step
  is stored in a
  \code{\link[BCRANK:BCRANKmatch-class]{BCRANKmatch-class}}
  object.  Several objects of this class are collected in a
  \code{\link[BCRANK:BCRANKresult-class]{BCRANKresult-class}} object}

\section{Objects from the Class}{
  Objects are not intended to be created directly but as a result from
  running \code{\link[BCRANK:bcrank]{bcrank}}.
}

\section{Slots}{
  \describe{
    \item{\code{searchPath}:}{a collection of \code{BCRANKmatch-class}
	  objects, ontaining all bcrank search steps from a start guess to a
	  locally optimal solution.}
	\item{\code{finalPWM}:}{position weight matrix containing all DNA
	  sites for the highest scoring consensus in this bcrank run.}
  }
}


\section{Methods}{
  \describe{
	\item{searchPath}{\code{signature(object = "BCRANKsearch", i=NULL)}:
	  If \code{i} is NULL, returns a data frame containing consensus and score
	  for the whole search path. Otherwise, the i'th
	  \code{BCRANKmatch-class} object in the search path is returned.}
    \item{pwm}{\code{signature(object = "BCRANKsearch", normalize=TRUE)}: Returns
	  the position weight matrix (pwm) for the highest scoring consensus in
	  this bcrank run. Matrix positions are between between 0 and 1 when \code{normalize}
	  is TRUE. When FALSE, the number of matching sequences is reported.}  
	  \item{plot}{\code{signature(x = "BCRANKsearch", y = "missing")}:
	 A plot method for the searchPath.}
 }
}

\author{Adam Ameur, \email{adam.ameur@lcb.uu.se}}

\seealso{
  \code{\link[BCRANK:bcrank]{bcrank}},
  \code{\link[BCRANK:BCRANKmatch-class]{BCRANKmatch-class}},
  \code{\link[BCRANK:BCRANKresult-class]{BCRANKresult-class}}
}

\keyword{classes}
