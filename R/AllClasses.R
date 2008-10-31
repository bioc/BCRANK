.onLoad <- function(lib, pkg) require("methods", quietly=TRUE)

## Class BCRANKmatch
setClass("BCRANKmatch",
         representation(## IUPAC consensus  
                        consensus="character",
                        ## Score
                        bcrankScore="numeric",
                        ## Vector with consensus occurences in ranked sequences.
                        ## For each region 0,1 or 2, where '2' represents 2 or more.
                        matchVec="numeric"))

## Class BCRANKsearch
setClass("BCRANKsearch",
         representation(
                        ## The resulting BCRANKmatch object
                        final="BCRANKmatch",
                        ## Position weight matrix for the final BCRANKmatch
                        finalPWM="matrix",
                        ## Total number of consensus occurences for the final BCRANKmatch
                        finalNrMatch="numeric",
                        ## A list of BCRANKmatch objects -- The whole search path 
                        searchPath="list",
                        ## Number of iterations from start guess to final
                        nrIterations="numeric"))

## Class BCRANKresult
setClass("BCRANKresult",
         representation(
                        ## The fasta file with ranked DNA sequences
                        fname = "character",
                        ## Number of sequences in fasta file
                        nrSeqs = "numeric",
                        ## The bcrank function call 
                        funCall="character",
                        ## A list of BCRANKsearch objects, ranked by their final scores
                        toplist="list",
                        ## A list of BCRANKsearch objects, ranked by their final scores
                        restarts="numeric"))
