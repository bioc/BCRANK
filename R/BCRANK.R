
###############################################################################
##
## File: BCRANK.R
##
## predicting binding site consensus from ranked DNA sequences
##
## Author: Adam Ameur, The Linnaeus Centre for Bioinformatics, Uppsala, Sweden
##
## Description: Implementation of an algorithm for detection of short DNA
##   sequences that are overrepresented in some part of the list. Starting from
##   some initial consensus DNA sequence coded in IUPAC symbols, the method uses
##   a heuristic search to improve the consensus until a local optimum is found.
## 
## Dependencies: Biostrings
##
################################################################################

####################### Help functions ######################### 

## IUPAC symbols for nucleotides
iupac <- c("N","A","C","G","T","R","Y","K","M","S","W","B","D","H","V")
names(iupac) <- c("ACGT","A","C","G","T","AG","CT","GT","AC","CG","AT","CGT","AGT","ACT","ACG")

## IUPAC symbols represented as a list
iupacAsList <- list("A"=c("A"),"C"=c("C"),"G"=c("G"),"T"=c("T"),"R"=c("A","G"),"Y"=c("C","T"),"K"=c("G","T"),"M"=c("A","C"),"S"=c("C","G"),"W"=c("A","T"),"B"=c("C","G","T"),"D"=c("A","G","T"),"H"=c("A","C","T"),"V"=c("A","C","G"),"N"=c("A","C","G","T"))

## Converts a string of IUPAC symbols to a consensus sequence
iupacToCons <- function(iupacString){

  iupacToBase <- function(base){
    if(base == "A" || base == "C" || base == "G" || base == "T")
      return(base)
    else{
      cons <- names(iupac[which(iupac == base)])
      return(paste("[",cons,"]",sep=""))
    }
  }

  iupacByPos <- strsplit(iupacString,"")[[1]]
  cons <- ""
  for(i in 1:length(iupacByPos)){
    cons <- paste(cons, iupacToBase(iupacByPos[i]),sep="")
  }
  return(cons)
} 


## Returns the reverse complement of a IUPAC or consensus sequence
revComp <- function(seq){
    
  baseComp <- array(NA,0)
  baseComp[c("A","C","G","T","R","Y","K","M","S","W","B","V","D","H","N","[","]","(",")")] <- c("T","G","C","A","Y","R","M","K","S","W","V","B","H","D","N","]","[",")","(")
  
  return(paste(rev(sapply(strsplit(toupper(seq),"")[[1]],function(x){
    if(!is.na(baseComp[x]))
      {
        return(baseComp[x])
      }else{
        return(x)
      }
  },USE.NAMES = FALSE)),collapse = ""))
}

## Trims leading and trailing white space from character strings.
trimWhiteSpace <- function(x){
    sub("[ \t\n\r]*$", "", sub("^[ \t\n\r]*", "", x))
}

## Removes all letters 'N' at the beginning and end of a string
trimNs <- function(consensus){
  return(sub("N*$", "", sub("^N*", "", consensus)))
}

## Adds n 'N' letters at the beginning and end of a string
addNs <- function(consensus,n){
  consensus <- paste(paste(rep("N",n),collapse=""),consensus,paste(rep("N",n),collapse=""),sep="")
  return(consensus)
}

## Reads sequences from a FASTA file and checks the file format..
seqFromFile <- function(fafile){

  ## Are only allowed to contain A,C,G,T 
  checkSequences <- function(sequences){
    sequenceMatch <- unlist(gregexpr("^[ACGT]*$",sequences))
    sequencesNotOK <- headers[which(sequenceMatch != 1)] 
    if(length(sequencesNotOK)>0){
      stop(paste("Error. Other letters than A,C,G,T in sequence:",sequencesNotOK[1],"'",sep=""))
    }
  }

  ## Read sequences using Biostings package
  res <- try(readFASTA(fafile), silent=TRUE)

  ## Report error
  if(class(res) == "try-error"){
    stop(paste("Error reading fasta file:",res))
  }
  
  seqs <- sapply(1:length(res), function(x){res[[x]]$seq})
  seqs <- trimWhiteSpace(seqs)

  headers <- lapply(1:length(res), function(x){res[[x]]$desc})
  
  names(seqs) <- headers

  ## Check for errors in DNA sequences
  checkSequences(seqs)
  
  return(seqs)
}


## Creates a position weight matrix (PWM) given sequences and a
## IUPAC consensus sequence. 
createPWM <- function(seqs, consensus, matchRevComp=TRUE){

  ## Create a position weight matrix as used by seqLogo
  createPWMhelp <- function(match, seqs, motifLength, isRevComp=FALSE){

    ## Detect all matching DNA sequences
    
    allTBS <- matrix(0, nrow=4,ncol=motifLength, dimnames=list(c("A","C","G","T"),1:motifLength))
    
    for(i in 1:length(match)){
      currentMatch <- match[[i]]
      TBS <- NULL

      if(currentMatch[1] > -1){
        seq <- seqs[i]
        for(j in 1:length(currentMatch)){
          posInRegion <- as.integer(currentMatch[j])
          matchLength <- attr(currentMatch,"match.length")[j]-1

          TBS <- as.character(substr(seq,posInRegion,posInRegion+matchLength))

          if(isRevComp==TRUE){
            TBS <- revComp(TBS)
          }
        }
        TBSasPWM <- sapply(strsplit(TBS,"")[[1]], function(x){
          if(x=="A"){
            return(c(1,0,0,0))
          }
          if(x=="C"){
            return(c(0,1,0,0))
          }
          if(x=="G"){
            return(c(0,0,1,0))
          }
          if(x=="T"){
            return(c(0,0,0,1))
          }
        })
        allTBS <- allTBS+TBSasPWM
      }
    }

    return(allTBS)
  }
    
  motifLength <- nchar(consensus)
    
  motifSequence <- iupacToCons(consensus)
   
  match <- gregexpr(motifSequence,seqs)
  pwm <- createPWMhelp(match, seqs, motifLength)
  
  if(matchRevComp == TRUE){
    RCmotifSequence <- revComp(motifSequence)
    matchRC <- gregexpr(RCmotifSequence,seqs)
    pwm_rc <- createPWMhelp(matchRC, seqs, motifLength,isRevComp=TRUE)
    pwm <- pwm+pwm_rc
  }

  return(pwm)

}


## Given a fasta file and a motif in IUPAC coding, the function makes a
## call to a C function that detects all occurences of the motif in the
## sequences. Number of sequences and the sequence lengths are required
## as arguments so that the C function can allocate enough memory.
matchConsensus <- function(seqs, motif, nrSeqs, seqLengths, maxOccurences=2){
  RCmotif <- revComp(motif)
  motifLength <- nchar(motif)
  RCmotif <- revComp(motif)
  
  hits <- .C("match_consensus", as.character(seqs), as.character(motif), as.character(RCmotif), as.integer(motifLength),  as.integer(nrSeqs), as.integer(seqLengths), as.integer(maxOccurences), result=as.integer(rep(0,nrSeqs)))$result
  
  return(hits)
}


## Given a fasta file and a motif in IUPAC coding, the function detects all
## occurences of the motif for all consensus sequences in the neighbourhood
## of the given motif.
hitsForNeighboursInC <- function(seqs, consensus, nrSeqs, seqLengths, silent){

  hitsForNeighboursByBase <- function(seqs, consensus, basePos, nrSeqs, seqLengths){

    ## Swap base at basePos to newBase
    swapBase<-function(consensus, basePos, newBase){
      consensusByPos <- strsplit(consensus,"")[[1]]
      
      prefix <- substr(consensus, 1, basePos-1)
      suffix <- substr(consensus, basePos+1, length(consensusByPos))
      
      return(paste(prefix,newBase,suffix,sep=""))
    }
    
    consensusByPos <- strsplit(consensus,"")[[1]]
    oldBase <-  substr(consensus, basePos, basePos)
    
    ## Consensus where base at position has been changed
    ## to A,C,G or T
    swappedA <- swapBase(consensus,basePos,"A")
    swappedC <- swapBase(consensus,basePos,"C")
    swappedG  <- swapBase(consensus,basePos,"G")
    swappedT <- swapBase(consensus,basePos,"T")
    
    basicHits <- list()
    basicHits[[swappedA]] <- matchConsensus(seqs, motif=swappedA, nrSeqs=nrSeqs, seqLengths=seqLengths)
    basicHits[[swappedC]] <- matchConsensus(seqs, motif=swappedC, nrSeqs=nrSeqs, seqLengths=seqLengths)
    basicHits[[swappedG]] <- matchConsensus(seqs, motif=swappedG, nrSeqs=nrSeqs, seqLengths=seqLengths)
    basicHits[[swappedT]] <- matchConsensus(seqs, motif=swappedT, nrSeqs=nrSeqs, seqLengths=seqLengths)
    
    hits <- list()
    
    newBases <- iupac[!(iupac %in% oldBase)]
    
    for(base in newBases){
      newCons <- swapBase(consensus, basePos, base)
      atomicBases <- iupacAsList[[base]]
      tmpHits <- rep(0,nrSeqs)
      for(baseTmp in atomicBases){
        consTmp <- swapBase(consensus, basePos, baseTmp)
        tmpHits <- tmpHits+basicHits[[consTmp]]
      }
      hits[[newCons]] <- tmpHits
    }
    
    return(hits)
  }

  consLength <- nchar(consensus)
  hits <- list()

  if(!silent){
    cat("Scanning sequences")
    flush.console()
  }
  
  for(basePos in 1:consLength){
    if(!silent){
      cat(".")
      flush.console()
    }
    hits <- append(hits,hitsForNeighboursByBase(seqs,consensus,basePos,nrSeqs,seqLengths))
  }

  if(!silent){
    cat("\n")
    flush.console()
  }
  
  return(hits)
}


#################### Neighbourhood #######################################

## Returns all consensus sequences in the neighbourhood of a given consensus
getBCRANKNeighbours <- function(consensus){
    
    consensusByPos <- strsplit(consensus,"")[[1]]
    
    neighbours <- array(NA)
    
    for(i in 1:length(consensusByPos)){
      prefix <- substr(consensus, 1, i-1)
      suffix <- substr(consensus, i+1, length(consensusByPos))
      
    newBases <- iupac[!(iupac %in% consensusByPos[i])]
      neighbours <- c(neighbours,paste(prefix,newBases,suffix,sep=""))
    }
  
    return(neighbours[which(!is.na(neighbours))])
  }

######################## Cost function ####################################

## Computes the BCRANK score of a consensus sequence
getBCRANKScore<-function(hitsVec, randOrder, consensus, use.P1, use.P2, nrPoints=25){
  
  ## Given the observed hits and the random orderings, the function returns
  ## the differences between the observed hits and each of the random
  ## orderings. Each difference is computed by comparing the cumulative
  ## hit functions at a number of uniformly distributed points in the interval.
  getDeviationsFromRandom <- function(hits, randOrder, nrPoints){

    points <- round(seq(1,length(hits),length.out=nrPoints))
    
    randHits <- sapply(1:ncol(randOrder), function(i){
      cumsum(hits[randOrder[,i]])[points]
    })

    obsHits <- cumsum(hits)[points]
    dists <- (obsHits - randHits)
        
    return(colSums(dists))
  }

  ## ###############
  ##
  ## Compute score
  ##
  ## ###############
  
  hitsLen <- length(hitsVec)
  hitsUnique <- rep(0,hitsLen)
  matchingIds <- which(hitsVec>0)
  nrHits <- length(matchingIds)

  ## Return negative score if motif is matching in none or all sequences
  if(nrHits == hitsLen || nrHits == 0){
    return(-100)
  }

  ## Make the hits vector conatin only 0's and 1's
  hitsUnique[matchingIds] <- 1
  hitsUnique <- hitsUnique/sum(hitsUnique)

  ## Compute significance score for the hits vector, by comparing to scores
  ## for random orderings using a t-test. 
  randDists <- getDeviationsFromRandom(hitsUnique, randOrder, nrPoints=nrPoints)
  tval <- abs(t.test(randDists,mu=0)$statistic)

  score <- tval
  
  ## #######################
  ##
  ## Compute score penalties.
  ##
  ## #######################
  
  ## Compute P1 - Penalty on non-specific bases

  if(use.P1){
    
    tmpConsensus <- trimNs(consensus)

    tmpConsensus <- gsub("A","X",tmpConsensus)
    tmpConsensus <- gsub("C","X",tmpConsensus)
    tmpConsensus <- gsub("G","X",tmpConsensus)
    tmpConsensus <- gsub("T","X",tmpConsensus)
    
    nrSpecific <- length(which(gregexpr("X",tmpConsensus)[[1]]!=-1))

    ## Avoid zero penalty for motifs with no A,C,G or T bases
    if(nrSpecific == 0){
      nrSpecific <- 0.5
    }

    P1 <-  nrSpecific/nchar(tmpConsensus)

    score <- score*P1

  }
    
  ## Compute P2 - Penalty on repetitive motifs
  if(use.P2){

    nrMultipleHits <- length(which(hitsVec>1))

    P2 <- 1-(nrMultipleHits/nrHits)

    score <- score*P2
  }
  
  return(score)
}


############################ Run search algorithm ###################################

## Runs the BCRANK search on a fasta file containing ranked DNA regions starting from
## an initial consensus sequence. The nrRandom parameter specifies the number of random
## re-orderings performed to calculate scores. Returns a BCRANKsearch object
bcrankRun <- function(seqs, start, nrRandom=500, silent=FALSE, makePlot=FALSE, do.search=TRUE, use.P1=TRUE, use.P2=TRUE){
    
  nrSeqs <- length(seqs)
  seqLengths <- as.numeric(unlist(lapply(seqs,nchar)))

  if(!silent){
    if(do.search){
      cat("\n**** Running BCRANK on",nrSeqs,"regions, starting from",start,"****\n")
    }
    else{
      cat("\n**** Computing BCRANK score for",start,"on",nrSeqs,"regions ****\n")
    }
    flush.console()
  }
  
  ## A list BCRANKmatch objects
  searchPath <- list()
    
  iteration <- 0 

  randOrder <-  sapply(1:nrRandom, function(x){
    sample(1:length(seqs))
  })

  ## A list with the consensus sequences that have alredy scores assigned to them.
  seenMotifs <- list()
  
  hits <- matchConsensus(seqs, motif=start, nrSeqs=nrSeqs, seqLengths=seqLengths)
  score <- getBCRANKScore(hits,randOrder,start, use.P1, use.P2)

  seenMotifs[[start]] <- score
    
  improvement <- TRUE
  iteration <- 0
  bestScore <- score
  bestCons <- start
  bestHits <- hits

  while(improvement){
    
    improvement <- FALSE
    iteration <- iteration+1
        
    searchPath[[iteration]] <- BCRANKmatch(as.character(bestCons),as.numeric(bestScore),as.numeric(bestHits))
    
    if(makePlot){
      ## Create a temporary BCRANKsearch object and plot it
      tmp_final <- searchPath[[iteration]]
      tmp_finalPWM <- matrix(NA)
      tmp_finalNrMatch <- 0
      tmp_searchObj <- BCRANKsearch(tmp_final, tmp_finalPWM, tmp_finalNrMatch, searchPath)
      plot(tmp_searchObj) 
    }

    if(do.search){
      if(!silent){
        cat("\nIteration ",iteration," - ",bestCons,": ",bestScore,"\n",sep="")
        flush.console()
      }
      
      ## Add flanking N's to the current best consensus
      bestCons <- addNs(bestCons,n=1)
      
      ## Generate neighbourhood around current best consensus
      neighbours <- getBCRANKNeighbours(bestCons)
      
      hitsNeigh <- hitsForNeighboursInC(seqs, bestCons, nrSeqs=nrSeqs, seqLengths=seqLengths, silent=silent)
      
      nrNeigh <- length(hitsNeigh)
      nrDone <- 0

      if(!silent){
        cat("Computing scores")
        flush.console()
      }
      
      for(consensus in names(hitsNeigh)){
        hits <- hitsNeigh[[consensus]]
        score <- NULL
        if(consensus %in% names(seenMotifs)){
          score <- seenMotifs[[consensus]]
        }
        else{
          score <- getBCRANKScore(hits,randOrder,consensus,use.P1,use.P2)
          seenMotifs[[consensus]] <- score
        }
        
        if(!is.nan(score)){
          if(score > bestScore){
            bestScore <- score
            bestCons <- trimNs(consensus)
            bestHits <- hits
            improvement <- TRUE
          }
        }
        
        nrDone <- nrDone+1
        
        if(!silent){
          ## Don't report all computed scores
          if((nrDone %% 10) == 0){
            cat(".")
            flush.console()
          }
        }
      }
      if(!silent){
        cat("\n")
        flush.console()
      }
    }
    ## Don't run the search
    else{
      improvement <- FALSE
    }
  }

  ## Create BCRANKmatch object for current best consensus and store in list
  final <- searchPath[[iteration]]
  finalPWM <- createPWM(seqs,consensus(final))
  finalNrMatch <- unique(colSums(finalPWM))
  
  ## Returns a BCRANKsearch object
  searchResult <- BCRANKsearch(final, finalPWM, finalNrMatch, searchPath)

  return(searchResult)
}


## Runs the BCRANK algorithm.
bcrank <- function(fafile, startguesses=c(), restarts=10, length=10, reorderings=500, silent=FALSE, plot.progress=FALSE, do.search=TRUE, use.P1=FALSE, use.P2=TRUE){
    
  ## A list containing BCRANKsearch objects
  BCRANKsearchResults <- list()
  
  ## Read sequences from file
  seqs <- seqFromFile(fafile)
    
  ## Generates a random consensus sequence in IUPAC encoding.
  getRandomInit <- function(length=10){
    return(paste(sample(iupac,length,replace=TRUE),collapse=""))
  }
  
  ## Multiple random search
  if(length(startguesses)==0){
    for(i in 1:restarts){
      BCRANKsearchResults[[i]] <- bcrankRun(seqs, nrRandom=reorderings, start=getRandomInit(length),silent=silent, makePlot=plot.progress, do.search=do.search, use.P1=use.P1, use.P2=use.P2)
    }
  }

  ## Search from given startguesses (coded in IUPAC)
  else{
    restarts <- length(startguesses)
    for(i in 1:restarts){
      BCRANKsearchResults[[i]] <- bcrankRun(seqs, nrRandom=reorderings, start=startguesses[i], silent=silent, makePlot=plot.progress, do.search=do.search, use.P1=use.P1, use.P2=use.P2)
    }
  }

  ## Return BCRANKresult object
  funCallString <- as.character(sys.call())
  result <- BCRANKresult(BCRANKsearchResults, length(seqs), fafile, funCallString, restarts)
  
  return(result)
}

