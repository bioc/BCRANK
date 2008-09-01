######################################################
##             INITIALIZE METHOD                    ##
######################################################

## Initialize method for class BCRANKsearch
setMethod("initialize", "BCRANKsearch",
          function(.Object,
                   final,
                   finalPWM,
                   finalNrMatch,
                   searchPath,
                   nrIterations){

            if(length(finalNrMatch)!=1){
              stop("finalNrMatch should have length 1")
            }
            
            callNextMethod()
          })

BCRANKsearch <- function(final, finalPWM, finalNrMatch, searchPath){

  nrIterations <- length(searchPath)
  
  new("BCRANKsearch", final=final, finalPWM=finalPWM, finalNrMatch=finalNrMatch, nrIterations=nrIterations, searchPath=searchPath)
}

######################################################
##             DISPLAY METHODS                      ##
######################################################

setMethod("show", "BCRANKsearch", function(object){
    
  cat("\nAn object of class \"BCRANKsearch\"\n\n")
  cat(paste("Search path, starting from ",consensus(searchPath(object,1)),":\n\n",sep=""))
  print(searchPath(object))
  cat("\nPosition weight matrix for search result (",consensus(final(object)),"):\n\n",sep="")
  print(pwm(object))
  cat("\n\nUse methods searchPath(object) and pwm(object) to access object slots.\n\n")
  
})

setMethod("plot", signature(x="BCRANKsearch",y="missing"), function(x,y,...){

  args <- list(...)

  ## Read arguments and use defaults if not set
  if(is.null(args[["cex.axis"]]))
    args[["cex.axis"]] <- 1.5
  if(is.null(args[["cex.lab"]]))
    args[["cex.lab"]] <- 1.5
  if (is.null(args[["xlab"]]))
    args[["xlab"]] <- "Ordered regions"
  if (is.null(args[["ylab"]]))
    args[["ylab"]] <- "Consensus occurrences"
  if (is.null(args[["type"]]))
    args[["type"]] <- "l"
  if (is.null(args[["lwd"]]))
    args[["lwd"]] <- 2
  if (is.null(args[["frame.plot"]]))
    args[["frame.plot"]] <- FALSE
  if (is.null(args[["col"]]))
    args[["col"]] <- rainbow(15)
  if (is.null(args[["plot.legend"]]))
    args[["plot.legend"]] <- TRUE 
  if (is.null(args[["scale"]]))
    args[["scale"]] <- FALSE

  scale <- FALSE
  if(args[["scale"]] == TRUE){
    scale <- TRUE
  }
  args[["scale"]] <- NULL
  
  plotLegend <- TRUE 
  if(args[["plot.legend"]] == FALSE){
    plotLegend <- FALSE
  }

  args[["plot.legend"]] <- NULL

  cols <- args[["col"]]
  
  nrIter <- nrIterations(x)
    
  sp <- searchPath(x)
  consensuses <- sp[,"Consensus"]
  bcrankScores <- sp[,"Score"]

  maxMatches <- max(sapply(1:nrIter,function(i){ nrMatch(searchPath(x,i)) }))
  
  for(i in 1:nrIter){
    currentMatch <- searchPath(x,i)
    currentHits <- matchVector(currentMatch)

    if(scale == TRUE){
      currentHits <- currentHits/sum(currentHits)
    }

    colIndex <- i

    if(i > length(cols)){
      colIndex <- (i %% length(cols))
    }
    
    currentCol <- cols[colIndex]
    
    if(i==1){
      if(scale == FALSE){
        if(is.null(args[["ylim"]])){
          args[["ylim"]] <- c(0,maxMatches)
        }
      }
      args[["x"]] <- cumsum(currentHits)
      args[["col"]] <- currentCol
      do.call("plot", args)
      args[["cex.axis"]] <- NULL
      args[["cex.lab"]] <- NULL
      args[["xlab"]] <- NULL
      args[["ylab"]] <- NULL
      args[["frame.plot"]] <- NULL
    }
    else{
      args[["x"]] <- cumsum(currentHits)
      args[["y"]] <- NULL
      args[["col"]] <- currentCol
      do.call("points",args)
    }
  }

  if(plotLegend == TRUE){
    legendPos <- NULL
    if(scale==TRUE){
      legendPos <- "bottomright"
    }
    else{
      legendPos <- "topleft"
    }
    legend(legendPos,cex=1,bty="n",legend=consensuses,col=cols,lty=2,lwd=5, xjust=0)
  }
  
})


######################################################
##             ACCESSOR METHODS                     ##
######################################################


setMethod("final", "BCRANKsearch", function(object){object@final})
setMethod("finalPWM", "BCRANKsearch", function(object){object@finalPWM})
setMethod("finalNrMatch", "BCRANKsearch", function(object){object@finalNrMatch})
setMethod("nrIterations", "BCRANKsearch", function(object){object@nrIterations})

setMethod("pwm", "BCRANKsearch", function(object, normalize=TRUE){
  pwm <- finalPWM(object)
  if(normalize){
    pwm <- pwm/finalNrMatch(object)
  }
  return(pwm)
})

setMethod("searchPath", "BCRANKsearch", function(object, i=NULL){
  BCRANKmatches <- object@searchPath

  ## Returns search path as a data.frame
  if(is.null(i)){
    pathConsensuses <- as.character(sapply(BCRANKmatches, function(x){ consensus(x) }))
    pathScores <- as.numeric(sapply(BCRANKmatches, function(x){ bcrankScore(x) }))

    return(data.frame("Consensus"=pathConsensuses,"Score"=pathScores))
  }
  ## Returns the i'th search element in the toplist
  else{
    if(i>0 && i<=nrIterations(object)){
      return(BCRANKmatches[[i]])
    }
    else{
      stop("Argument i is out of range.")
    }
  }
})
          
