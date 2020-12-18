genemodel = function(sym, genome="hg19") {
 stopifnot(genome=="hg19")
 require(TxDb.Hsapiens.UCSC.hg19.knownGene)
 txdb = TxDb.Hsapiens.UCSC.hg19.knownGene
 require(org.Hs.eg.db)
 num = AnnotationDbi::get(sym, revmap(org.Hs.egSYMBOL))
 exonsBy(txdb, by="gene")[[num]]
}


setClass("bindingScores",
 representation(sco = "numeric", tfn="character"))

makebs = function(tfn) {
 require(harbChIP)
 if (!exists("harbChIP")) data(harbChIP)
 tmp = exprs(harbChIP)[,tfn]
 cl = na.omit(tmp)
 ans = as.numeric(cl)
 names(ans) = names(cl)
 new("bindingScores", sco=sort(ans, decreasing=TRUE), tfn=tfn)
}

setMethod("show", "bindingScores", function(object){
 cat("bindingScores for ", object@tfn, "\ntop 5:\n")
 print( head( object@sco ))
})

setGeneric("QQnorm", function(x, tx=force, ...) standardGeneric("QQnorm"))
setMethod("QQnorm", "bindingScores", function(x, ...) {
  sco = x@sco
#  require(parody)
  out = calout.detect(tx(sco), ...)
  lowout = tx(sco[ max(out$ind[out$ind < length(sco)/2]) ])
  qstr = qqnorm( tx(sco), 
     main=paste0("QQnormal: binding scores for ", x@tfn))
  abline(h=lowout)
  invisible(out)
})

setGeneric("boundGenes", function(x, ...)
  standardGeneric("boundGenes"))
setMethod("boundGenes", "bindingScores", function(x, tx=force, ...) {
#  require(parody)
  out = calout.detect(tx(x@sco), ...)
  names(x@sco)[out$ind]
})

setMethod("[", "bindingScores", function(x, i, j, ..., drop=FALSE) {
  x@sco = x@sco[i]
  x
})
    
