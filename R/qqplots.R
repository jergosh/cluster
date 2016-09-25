ggd.qqplot = function(pvector, main=NULL, ...) {
  o = -log10(sort(pvector,decreasing=F))
  e = -log10( 1:length(o)/length(o) )
  
  plot(e,o,pch=19,cex=1, main=main, ...,
       xlim=c(0,max(e)), ylim=c(0,max(o)),
       ann=FALSE)
  mtext(expression(Expected~~-log[10](italic(p))), side=1, line=2.5, cex=0.7)
  mtext(expression(Observed~~-log[10](italic(p))), side=2, line=1.75, cex=0.7)
  
  lines(e,e,col="red")
}

ggd.qqplot.mult = function(pvector, factor, cols, all=F, main=NULL, ...) {
  es <- list()
  os <- list()
  max_e <- 0.0
  max_o <- 0.0
  for (l in levels(factor)) {
    o <- -log10(sort(pvector[factor == l], decreasing=F))
    e <- -log10(1:length(o)/length(o))
    
    os[[l]] <- o
    es[[l]] <- e
    
    max_e <- max(c(max_e, e))
    max_o <- max(c(max_o, o))                
  }
  
  
  plot(NA, main=main, ...,
       xlim=c(0, max_e), ylim=c(0, max_o),
       ann=FALSE)
  mtext(expression(Expected~~-log[10](italic(p))), side=1, line=2.5)
  mtext(expression(Observed~~-log[10](italic(p))), side=2, line=2.0)
  
  
  if (all) {
    o <- -log10(sort(pvector, decreasing=F))
    e <- -log10(1:length(o)/length(o))
    
    points(e, o, pch=1, cex=1, col="black") 
  }
  
  for (l in levels(factor)) {
    points(es[[l]], os[[l]], pch=1, cex=1, col=cols[which(levels(factor) %in% l)])
  }
  
  lines(c(0, max_e), c(0, max_e), col="black")
}
