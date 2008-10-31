######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class BCRANKmatch
setMethod("initialize", "BCRANKmatch",
          function(.Object,
                   consensus,
                   bcrankScore,
                   matchVec){
            
            ## Check that consensus is in IUPAC coding
            iupacMatch <- gregexpr("^[ACGTRYKMSWBVDHN]*$",consensus)[[1]]
            
            if(attr(iupacMatch,"match.length") < 0){
              stop(paste(consensus,"not in IUPAC coding."))
            }
            
            ##callNextMethod()

            .Object@consensus <- consensus
            .Object@bcrankScore <- bcrankScore
            .Object@matchVec <- matchVec
            .Object
            
          })

BCRANKmatch <- function(consensus, bcrankScore, matchVec){
  new("BCRANKmatch", consensus=consensus, bcrankScore=bcrankScore, matchVec=matchVec)
}

######################################################
##             DISPLAY METHODS                      ##
######################################################

setMethod("show", "BCRANKmatch", function(object){
    
  cat("\nAn object of class \"BCRANKmatch\"\n\n")
  cat("Consensus:",consensus(object),", bcrankScore:",bcrankScore(object))
  
  cat("\n\nUse methods consensus(object), bcrankScore(object) and matchVector(object) to access object slots.\n\n")
  
})


######################################################
##             ACCESSOR METHODS                     ##
######################################################

setMethod("consensus", "BCRANKmatch", function(object) object@consensus)
setMethod("bcrankScore", "BCRANKmatch", function(object) object@bcrankScore)

setMethod("matchVector", "BCRANKmatch", function(object){
  matchVec <- object@matchVec
  matchVec[which(matchVec>0)] <- 1
  return(matchVec)
})

setMethod("nrMatch", "BCRANKmatch", function(object){
  return(length(which(matchVector(object) != 0)))
})
