library('minpack.lm')
#source('~/Desktop/Rhythmicity Code and Output/fitSinCurve.R')
expr.c<-read.csv("~/Desktop/Rhythmicity Code and Output/Example_data.csv", row.names = 1)
clin_ctrl<-read.csv('~/Desktop/Rhythmicity Code and Output/Example_clinical.csv', row.names = 1)
n<-nrow(expr.c)
Symbols<-row.names(expr.c)
observed_para_c <- data.frame(genes=Symbols,A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
tod_c<-clin_ctrl$TOD
for (i in 1:n) {
  out <- fitSinCurve(xx=tod_c,observed=as.numeric(expr.c[i,]))
  observed_para_c[i,-1] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  if(i%%1000==0) print(i)
}
observed_para_c_sorted<-observed_para_c[order(observed_para_c$R2, decreasing = TRUE),]

#Generate Scatter plots 
#setwd("~/Desktop/Rhythmicity Code and Output/")
top.control<-observed_para_c[order(observed_para_c$R2, decreasing = TRUE)[1:4],]
index<-as.numeric(row.names(top.control))
labels1<-clin_ctrl$Site
specInfo<-"Control"
for(i in 1:nrow(top.control)){
  agene <- as.character(top.control$genes[i])
  fileName <- paste0('Top_Control',agene,'.pdf')
  pdf(fileName)
  circadianDrawing(tod=tod_c, expr=expr.c[index[i],], apar=top.control[i,],labels=labels1, specInfo=specInfo)
  dev.off()
}

##Generate p/q values 
#setwd("~/Desktop/Rhythmicity Code and Output/")
#system("mkdir -p nullFolder")
setwd("nullFolder")
groupName <- 'control'
thisData <- expr.c
library(doParallel)
registerDoParallel()
B<-10 #we use 1000 permutations in the paper, but to save time and computation power, 
#the results in this folder reflect 10 permutations
result <- foreach(b = 1:B) %dopar% {
  print(b)	
  library(minpack.lm)
  source('~/Desktop/Rhythmicity Code and Output/fitSinCurve.R')
  
  null_pare <- data.frame(A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
  null_para_file <- paste('null_',groupName,'_',b,'.rdata',sep='')
  
  set.seed(b)
  shuffleTOD <- sample(tod_c)
  
  for (i in 1:n) {
    out <- fitSinCurve(xx=shuffleTOD,observed=unlist(thisData[i,]))
    null_pare[i,] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  }		
  save(null_pare,file=null_para_file)	
}


null_pare_A <- matrix(0,n,B)
null_pare_phase <- matrix(0,n,B)
null_pare_offset <- matrix(0,n,B)
null_pare_peak <- matrix(0,n,B)
null_pare_R2 <- matrix(0,n,B)

for(b in 1:B){
  print(b)
  file11 <- paste('null_',groupName,'_',b,'.rdata',sep='')
  
  load(file11)
  null_pare_A[,b] <- null_pare$A
  null_pare_phase[,b] <- null_pare$phase
  null_pare_offset[,b] <- null_pare$offset
  null_pare_peak[,b] <- null_pare$peak
  null_pare_R2[,b] <- null_pare$R2		
}

null_para <- list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
                  null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)

null_para_file <- paste('null_',groupName,'.rdata',sep='')
save(null_para,file=null_para_file)

null_para <- get(load("~/Desktop/Rhythmicity Code and Output/nullFolder/null_control.rdata"))
para_R2_pool <- c(observed_para_c$R2,null_para$null_para_R2)
R2Rank_para <- 1 - (rank(para_R2_pool)[1:length(observed_para_c$R2)] - 0.5)/length(para_R2_pool)
observed_para_c$pvalue <- R2Rank_para
observed_para_c$qvalue <- p.adjust(observed_para_c$pvalue, 'BH')
observe_para_c_sorted<-observed_para_c[order(observed_para_c$pvalue),]
#setwd("~/Desktop/Rhythmicity Code and Output/")
write.csv(observe_para_c_sorted, "Example_result.csv")

##Rhythmicity gain/loss example 
#Split data into two groups: For example we will split the controls into younger and older 
old_index<-which(clin_ctrl$Age>=50)
young_index<-which(clin_ctrl$Age<50)

expr.old<-expr.c[,old_index]
expr.young<-expr.c[,young_index]

tod_o<-clin_ctrl$TOD[old_index]
tod_y<-clin_ctrl$TOD[young_index]

observed_para_y <- data.frame(genes=Symbols,A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))

for (i in 1:n) {
  out <- fitSinCurve(xx=tod_y,observed=as.numeric(expr.young[i,]))
  observed_para_y[i,-1] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  if(i%%1000==0) print(i)
}

observed_para_o <- data.frame(genes=Symbols,A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))

for (i in 1:n) {
  out <- fitSinCurve(xx=tod_o,observed=as.numeric(expr.old[i,]))
  observed_para_o[i,-1] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  if(i%%1000==0) print(i)
}

#Generate Null Distributions: 
#setwd('~/Desktop/Rhythmicity Code and Output/nullFolder')
library(doParallel)
groupName <- 'old'
thisData <- expr.old
B<-10
result <- foreach(b = 1:B) %dopar% {
  print(b)	
  library(minpack.lm)
  source('~/Desktop/Rhythmicity Code and Output/fitSinCurve.R')
  
  null_pare <- data.frame(A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
  null_para_file <- paste('null_',groupName,'_',b,'.rdata',sep='')
  
  set.seed(b)
  shuffleTOD <- sample(tod_o)
  
  for (i in 1:n) {
    out <- fitSinCurve(xx=shuffleTOD,observed=unlist(thisData[i,]))
    null_pare[i,] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  }		
  save(null_pare,file=null_para_file)	
}

null_pare_A <- matrix(0,n,B)
null_pare_phase <- matrix(0,n,B)
null_pare_offset <- matrix(0,n,B)
null_pare_peak <- matrix(0,n,B)
null_pare_R2 <- matrix(0,n,B)

for(b in 1:B){
  print(b)
  file11 <- paste('null_',groupName,'_',b,'.rdata',sep='')
  
  load(file11)
  null_pare_A[,b] <- null_pare$A
  null_pare_phase[,b] <- null_pare$phase
  null_pare_offset[,b] <- null_pare$offset
  null_pare_peak[,b] <- null_pare$peak
  null_pare_R2[,b] <- null_pare$R2		
}

null_para <- list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
                  null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)

null_para_file <- paste('null_',groupName,'.rdata',sep='')
save(null_para,file=null_para_file)

###REPEAT NULL GENERATION FOR Young###
#setwd('~/Desktop/Rhythmicity Code and Output/nullFolder')
library(doParallel)
groupName <- 'young'
thisData <- expr.young
B<-10
result <- foreach(b = 1:B) %dopar% {
  print(b)	
  library(minpack.lm)
  source('~/Desktop/Rhythmicity Code and Output/fitSinCurve.R')
  
  null_pare <- data.frame(A=numeric(n), phase=numeric(n), offset=numeric(n), peak=numeric(n), R2=numeric(n))
  null_para_file <- paste('null_',groupName,'_',b,'.rdata',sep='')
  
  set.seed(b)
  shuffleTOD <- sample(tod_y)
  
  for (i in 1:n) {
    out <- fitSinCurve(xx=shuffleTOD,observed=unlist(thisData[i,]))
    null_pare[i,] <- c(out$A, out$phase, out$offset, out$peak, out$R2)
  }		
  save(null_pare,file=null_para_file)	
}

null_pare_A <- matrix(0,n,B)
null_pare_phase <- matrix(0,n,B)
null_pare_offset <- matrix(0,n,B)
null_pare_peak <- matrix(0,n,B)
null_pare_R2 <- matrix(0,n,B)

for(b in 1:B){
  print(b)
  file11 <- paste('null_',groupName,'_',b,'.rdata',sep='')
  
  load(file11)
  null_pare_A[,b] <- null_pare$A
  null_pare_phase[,b] <- null_pare$phase
  null_pare_offset[,b] <- null_pare$offset
  null_pare_peak[,b] <- null_pare$peak
  null_pare_R2[,b] <- null_pare$R2		
}

null_para <- list(null_para_A=null_pare_A, null_para_phase=null_pare_phase, null_para_offset=null_pare_offset, 
                  null_para_peak=null_pare_peak, null_para_R2=null_pare_R2)

null_para_file <- paste('null_',groupName,'.rdata',sep='')
save(null_para,file=null_para_file)

R2change <- observed_para_o$R2 - observed_para_y$R2 
Ashift <- abs(abs(observed_para_o$A) - abs(observed_para_y$A))
phaseshift0 <- observed_para_o$phase - observed_para_y$phase
phaseshift <- pmin(abs(phaseshift0), 24 - abs(phaseshift0))
intshift <- abs(observed_para_o$offset - observed_para_y$offset)

R2changeNULL <- R2change
AshiftNULL <- Ashift
phaseshiftNULL <- phaseshift
intshiftNULL <- intshift

B <- 10
for(b in 1:B){
  if(b%%50 == 0) print(b)
  file_old <- paste0("~/Desktop/Rhythmicity Code and Output/nullFolder/",'null_old_',b,'.rdata')
  file_young <- paste0("~/Desktop/Rhythmicity Code and Output/nullFolder/",'null_young_',b,'.rdata')
  null_para_old <- get(load(file_old))
  null_para_young <- get(load(file_young))
  
  aR2change <- null_para_old$R2 - null_para_young$R2
  aAshift <- abs(abs(null_para_old$A) - abs(null_para_young$A))
  aphaseshift <- pmin(abs(null_para_old$phase - null_para_young$phase),24 - abs(null_para_old$phase - null_para_young$phase))
  aintshift <- abs(null_para_old$offset - null_para_young$offset)
  
  R2changeNULL <- c(R2changeNULL,aR2change)
  AshiftNULL <- c(AshiftNULL,aAshift)
  phaseshiftNULL <- c(phaseshiftNULL,aphaseshift)
  intshiftNULL <- c(intshiftNULL,aintshift)
}

p <- nrow(observed_para_o)

R2gainPvalue <- 1-rank(R2changeNULL)[1:p]/length(R2changeNULL) 
R2losePvalue <- rank(R2changeNULL)[1:p]/length(R2changeNULL) 
AshiftPvalue <- 1 - rank(AshiftNULL)[1:p]/length(AshiftNULL)
phasePvalue <- 1 - rank(phaseshiftNULL)[1:p]/length(phaseshiftNULL)
intshiftPvalue <- 1 - rank(intshiftNULL)[1:p]/length(intshiftNULL)

result2<-data.frame(cbind(R2gainPvalue, R2losePvalue, AshiftPvalue, phasePvalue, intshiftPvalue))
row.names(result2)<-observed_para_o$genes
result2_sorted<-result2[order(result2$R2gainPvalue, decreasing = FALSE), ]
#setwd("~/Desktop/Rhythmicity Code and Output/")
write.csv(result2_sorted, "Example_Result2.csv")
