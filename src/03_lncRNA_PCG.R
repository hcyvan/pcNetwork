library(dplyr)
source('./lib/globals.R')
source('./lib/helpers.R')
load('./cache/diff.qlf.2877.cor.pvalue.rda')

load('./cache/diff.qlf.2877.pairs.rda') # diff.cor.pairs
diff.gene <- helper.get.lncRNA.PCG()

lncRNA <- filter(diff.gene, GeneType%in%config$lncRNA)$GeneID
PCG <- filter(diff.gene, GeneType%in%config$PCGs)$GeneID

lncRNA2PCG <- list()
i <- 1
for (lr in lncRNA) {
  cat(paste(i, lr, '\n'))
  i <- i+1
  tmp <- list()
  data <- filter(diff.cor.pairs, Var1==lr & type2%in%config$PCGs)
  tmp[['data']] <- data
  tmp[['pcg']] <- as.vector(filter(data, significant)$Var2)
  tmp[['pcg.r']] <- length(tmp[['pcg']] )
  
  data.0.7 <- filter(data, significant&abs(r)>=0.7)
  tmp[['pcg.0.7']] <- as.vector(data.0.7$Var2)
  tmp[['pcg.0.7.r']] <- as.vector(data.0.7$r)
  
  data.0.8 <- filter(data, significant&abs(r)>=0.8)
  tmp[['pcg.0.8']] <- as.vector(data.0.8$Var2)
  tmp[['pcg.0.8.r']] <- as.vector(data.0.8$r)
  
  lncRNA2PCG[[lr]] <- tmp
}

save(lncRNA2PCG, file='./cache/lncRNA2PCG.rda')
load('./cache/lncRNA2PCG.rda')

lapply(lncRNA2PCG, function(rna){
  rna$pcg.0.7
}) -> pcgs
sapply(pcgs, length)
# pcg
#
#


Reduce(union,pcgs) ->all

lapply(lncRNA2PCG, function(rna){
  as.vector(rna$data$r[match(all,rna$data$Var2)])
}) -> rna.pcg.r
mat <- matrix(unlist(rna.pcg.r), ncol = length(rna.pcg.r))
mat <- mat[,sapply(pcgs, length)!=0]

heatmap(mat)
