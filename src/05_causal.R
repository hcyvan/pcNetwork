
library(ctest)
f <- function(x,y,p=100){
  out <- try(c(unlist(permANM(x,y,p)[c(1,4,5)])))
  if(!is.numeric(out)){
    return(c(NA,0,NA))
  }
  return(out)
}
raw_pairs <- select(filter(diff.cor.pairs,significant&abs(r)>=.8),Var1,Var2)
# pairs <- raw_pairs[1:100,]
pairs <- raw_pairs; dim(pairs)

###################################
# First Round
###################################

fdata <- fpkm.data
system.time(out <- sapply(1:nrow(pairs),function(i){
  if(i%%100==1){print(paste(Sys.time(),i))}
  xi1 <- fdata[fdata$GeneID==pairs[i,1],-1]
  xi2 <- fdata[fdata$GeneID==pairs[i,2],-1]
  f(as.numeric(xi1),as.numeric(xi2),1)
}))
sum(out[2,]==0)
outbk <- out
save(outbk,file='test_out.rda')
#0 no direction, 1;-1 directed

###################################
# Rolling
###################################

pairs <- pairs[out[2,]==0,]
system.time(out <- sapply(1:nrow(pairs),function(i){
  if(i%%100==1){print(paste(Sys.time(),i))}
  xi1 <- fdata[fdata$GeneID==pairs[i,1],-1]
  xi2 <- fdata[fdata$GeneID==pairs[i,2],-1]
  f(as.numeric(xi1),as.numeric(xi2),10)
}))
sum(out[2,]==0)

pairs <- pairs[out[2,]==0,]
system.time(out2 <- sapply(1:nrow(pairs),function(i){
  if(i%%100==1){print(paste(Sys.time(),i))}
  xi1 <- fdata[fdata$GeneID==pairs[i,1],-1]
  xi2 <- fdata[fdata$GeneID==pairs[i,2],-1]
  f(as.numeric(xi1),as.numeric(xi2),100)
}))
sum(out2[2,]==0)

###################################
# 100 pairs
###################################
p100 <- data.frame(pairs,t(out2))
p100 <- filter(p100,P_no_causation==0&!is.na(dir))
key.p100 <- paste(p100$Var1,p100$Var2)
key.dcp <- paste(diff.cor.pairs$Var1,diff.cor.pairs$Var2)
p100 <- data.frame(p100,diff.cor.pairs[match(key.p100,key.dcp),-1:-2])

