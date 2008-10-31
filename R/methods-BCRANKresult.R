######################################################
##             INITIALIZE METHOD                    ##
######################################################

BCRANKresult <- function(BCRANKsearchList,nrSeqs,fname, funCall, restarts){

  scores <- sapply(1:length(BCRANKsearchList), function(x){bcrankScore(final(BCRANKsearchList[[x]]))})
  
  ordering <- order(scores,decreasing=TRUE)

  toplistOrdered <- lapply(ordering, function(x){BCRANKsearchList[[x]]})
  
  new("BCRANKresult", toplist=toplistOrdered, nrSeqs=nrSeqs, fname=fname, funCall=funCall, restarts=restarts)
}

######################################################
##             DISPLAY METHODS                      ##
######################################################

setMethod("show", "BCRANKresult", function(object){
    
  cat("\nAn object of class \"BCRANKresult\"\n\n")
  cat(paste("Top",restarts(object),"DNA motifs predicted by BCRANK:\n\n"))
  
  print(toptable(object))

  cat("\n\nUse methods toptable(object) and fname(object) to access object slots.\n\n")
  
})

######################################################
##             ACCESSOR METHODS                     ##
######################################################


setMethod("restarts", "BCRANKresult", function(object){object@restarts})
setMethod("toplist", "BCRANKresult", function(object){object@toplist})
setMethod("fname", "BCRANKresult", function(object){object@fname})
setMethod("funCall", "BCRANKresult", function(object){object@funCall})
setMethod("nrSeqs", "BCRANKresult", function(object){object@nrSeqs})

setMethod("toptable", "BCRANKresult", function(object, i=NULL){
  BCRANKsearches <- toplist(object)

  ## Returns toplist as a data.frame
  if(is.null(i)){
    topConsensuses <- as.character(sapply(BCRANKsearches, function(x){ consensus(final(x)) }))
    topScores <- as.numeric(sapply(BCRANKsearches, function(x){ bcrankScore(final(x)) }))

    return(data.frame("Consensus"=topConsensuses,"Score"=topScores))
  }
  ## Returns the i'th search element in the toplist
  else{
    if(i>0 && i<=restarts(object)){
      return(BCRANKsearches[[i]])
    }
    else{
      stop("Argument i is out of range.")
    }
  }
                        
})
