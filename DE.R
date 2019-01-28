##DE Analysis 
expr.c<-read.csv("~/Desktop/Rhythmicity Code and Output/Example_data.csv", row.names = 1)
clin_ctrl<-read.csv('~/Desktop/Rhythmicity Code and Output/Example_clinical.csv', row.names = 1)
n<-nrow(expr.c)
Symbols<-row.names(expr.c)

TOD<-clin_ctrl$TOD
night.index.neg<-which(TOD < 0)
night.index.pos<-which(TOD<=18 & TOD>12)
night.index<-c(night.index.neg,night.index.pos) #32
morn.index<-which(TOD>=0 & TOD<=12) #60
tod_index<-rep(1,104)
tod_index<-replace(tod_index, morn.index, 0)
clin_ctrl$tod_index<-tod_index

covariates0<-c("Gender","PMI", "Age", "Ethnicity", "Site")
covariables3 <- clin_ctrl[,covariates0]
UniqueID <-clin_ctrl$Site
VariableListOne = combn(colnames(covariables3),1)
VariableListTwo = combn(colnames(covariables3),2)
Diagnosis<-as.factor(tod_index)
new.data<-expr.c

library(snowfall) #parallel computing to speed up process 
set.seed(15213)
L3startTime = Sys.time()
sfInit(parallel=TRUE,type="SOCK", cpus=24) 
GeneIndeces = 1:nrow(new.data)
sfExport("new.data") 
sfExport("covariables3") 
sfExport("Diagnosis") 
sfExport("UniqueID") 
sfExport("VariableListOne") 
sfExport("VariableListTwo") 
sfExport("GeneIndeces") 

indRIM<-c()
sfExport("bestModelSelection") ##Very long, just running first 100 genes 
indRIM<-sfLapply(1:100,function(x) try(bestModelSelection(x,new.data,covariables3,Diagnosis,UniqueID,VariableListOne,VariableListTwo),silent=TRUE ) ) 
sfStop()

#Permutation to correct for p-value bias 
rawdata8<-new.data[1:100,]
B<-10
#Corrected 
#setwd("~/Desktop/DE Code and Output/")
#system("mkdir NullData")
setwd("NullData")
for(b in 1:B){
  print(b)
  set.seed(b)
  rawdata8_b <- rawdata8[,sample(ncol(rawdata8))]
  result_PVL3_SZ <- replicate(nrow(rawdata8_b),list())
  names(result_PVL3_SZ) <- rownames(rawdata8_b)
  sfExport("rawdata8_b")
  result_PVL3_SZ<-sfLapply(1:nrow(rawdata8_b),function(x) try(bestModelSelection(x,rawdata8_b,covariables3,Diagnosis,UniqueID,VariableListOne,VariableListTwo),silent=TRUE ) ) 
  afile <- paste('DE_null_',b,'.rdata',sep='')
  save(result_PVL3_SZ,file=afile)
}

RIM.NULL<-c()
for(b in 1:10){
  print(b)
  aNull<-get(load(paste("DE_null_", b,".rdata", sep = "")))
  RIM.NULL[[b]]<-aNull
}

unlistRIM<-unlist(RIM.NULL)
names.pval<-grep("lrt.pvalue", names(unlistRIM))
pval.null<-unname(as.numeric(unlistRIM[names.pval]))
unlistRIM2<-unlist(indRIM)
names.pval<-grep("lrt.pvalue", names(unlistRIM2))
pval.biased<-unname(as.numeric(unlistRIM2[names.pval]))

p.corr<-c()
numerator<-c()
denominator<-length(pval.null)
for(i in 1:length(pval.biased)){
  numerator[i]<-length(which(pval.null<pval.biased[i])) #empirical p valuse calculation
  p.corr[i]<-numerator[i]/denominator
  print(i) 
} 
bh.q<-p.adjust(p.corr, "BH")
library(qvalue)
q<-qvalue(p.corr)$qvalues
result3<-data.frame(cbind(pval.biased, p.corr, bh.q, q))
row.names(result3)<-row.names(new.data[1:100,])
#setwd("~/Desktop/DE Code and Output/")
#write.csv(result3, "Example_result3.csv")
