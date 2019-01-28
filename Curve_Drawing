#Circadian Drawing
circadianDrawing <- function(tod, expr, apar, labels, specInfo=NULL){	
  getPred <- function(parS, xx) {	
    parS$A * sin(2*pi/24 * (xx + parS$phase)) + parS$offset
  }
  
  geneName <- apar$genes
  peak <- round(apar$peak)
  if(peak==18) peak <- -6
  #pvalue <- signif(apar$pvalue,3)
  
  # amain <- paste('PV L3 healthy\n',geneName,':',probeName,'\n','p-value =',apvalue,sep='')
  amain <- paste(specInfo,', ',geneName,': ','; peak = ',peak,sep='')
  
  times <- seq(-6,18,0.1)
  pred <- getPred(apar,times)
  
  labelColor <- as.numeric(factor(labels))
  
  plot(tod,expr,col=labelColor, pch=16,cex=2,
       main=amain,xlim=c(-6,18),
       xlab='TOD',ylab='Expression')
  smoothingSpline = smooth.spline(times, pred, spar=0.35)
  lines(smoothingSpline,col='red',lwd=4)
  box(which = "plot", lty = "solid",lwd=3)	
  #legend('topright',legend=unique(labels),col=unique(labelColor),pch=16,cex=2)
}
