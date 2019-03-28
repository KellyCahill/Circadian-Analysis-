## need package minpack.lm
fitSinCurve <- function(xx, observed, parStart = list(A=3,phase=0, offset=0)){
	getPred <- function(parS, xx) {	
		parS$A * sin(2*pi/24 * (xx + parS$phase)) + parS$offset
	}

	residFun <- function(p, observed, xx) observed - getPred(p,xx)

	nls.out <- nls.lm(par=parStart, fn = residFun, observed = observed,	xx = xx)
	
	apar <- nls.out$par
	
	A0 <- apar$A
	asign <- sign(A0)
	## restrict A > 0
	A <- A0 * asign
	#phase <- (round(apar$phase) + ifelse(asign==1,0,12)) %% 24 
	phase <- (apar$phase + ifelse(asign==1,0,12)) %% 24 
	offset <- apar$offset
	
	peak <- (12 * sign(A0) - 6 - phase) %%24
	if(peak > 18) peak = peak - 24
	
	SSE <- sum(nls.out$fvec^2)
	SST <- sum((observed - mean(observed))^2)
	R2 <- 1 - SSE/SST
	res <- list(A=A, phase=phase, offset=offset, peak=peak, R2=R2)
	res
}
