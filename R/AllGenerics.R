## GENERIC ACCESSOR METHODS:

## BCRANKmatch
setGeneric("consensus", function(object) standardGeneric("consensus"))
setGeneric("bcrankScore", function(object) standardGeneric("bcrankScore"))
setGeneric("matchVector", function(object) standardGeneric("matchVector"))
setGeneric("nrMatch", function(object) standardGeneric("nrMatch"))

## BCRANKsearch
setGeneric("final", function(object) standardGeneric("final"))
setGeneric("finalPWM", function(object) standardGeneric("finalPWM"))
setGeneric("finalNrMatch", function(object) standardGeneric("finalNrMatch"))
setGeneric("nrIterations", function(object) standardGeneric("nrIterations"))
setGeneric("searchPath", function(object, i=NULL) standardGeneric("searchPath"))
setGeneric("pwm", function(object, normalize=TRUE) standardGeneric("pwm"))

## BCRANKresult
setGeneric("fname", function(object) standardGeneric("fname"))
setGeneric("nrSeqs", function(object) standardGeneric("nrSeqs"))
setGeneric("funCall", function(object) standardGeneric("funCall"))
setGeneric("toplist", function(object) standardGeneric("toplist"))
setGeneric("restarts", function(object) standardGeneric("restarts"))
setGeneric("toptable", function(object, i=NULL) standardGeneric("toptable"))
