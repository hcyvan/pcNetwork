library(dplyr)
source('./lib/globals.R')
source('./lib/helpers.R')

############################################### lncRNA-pcg
load('./cache/diff.qlf.2877.pairs.rda') # diff.cor.pairs

get.lncRNA.2.PCG <- function(s=0) {
  data <- filter(diff.cor.pairs, type1%in%config$lncRNA & type2%in%config$PCGs&significant&abs(r)>=s)
  lncRNA.2.PCG <- lapply(split(data,as.vector(data$Var1)), function(x){as.vector(x$Var2)})
  lncRNA.2.PCG
}

lncRNA.2.PCG.plot <- function(s) {
  sapply(get.lncRNA.2.PCG(s), function(x){length(x)})%>%sort->tmp
  tmp%>%barplot(main=s)
  summary(tmp)
}
par(mfrow=c(3,3))
lncRNA.2.PCG.plot(0.1)
lncRNA.2.PCG.plot(0.2)
lncRNA.2.PCG.plot(0.3)
lncRNA.2.PCG.plot(0.4)
lncRNA.2.PCG.plot(0.5)
lncRNA.2.PCG.plot(0.6)
lncRNA.2.PCG.plot(0.7)
lncRNA.2.PCG.plot(0.8)
lncRNA.2.PCG.plot(0.9)

################################################# tf-pcg
get.tf.2.PCG.from.enricher <- function() {
  tf.enricher <- read.delim('./data/enricher.all.bg0.new.csv', sep = ',', stringsAsFactors = F)
  tf.enricher.list <- setNames(split(tf.enricher, seq(nrow(tf.enricher))), tf.enricher$detail.ID)
  
  load('./cache/biomart.symbol.biotype.rda')
  
  lapply(tf.enricher.list, function(tf){
    symbol <- str_split(tf$detail.geneID, '/')[[1]]
    filter(biomart.symbol.biotype, hgnc_symbol%in%symbol)$ensembl_gene_id
  })
}

get.tf.2.PCG.from.fimo <- function() {
  fimo <- read.delim('./data/enricher/fimo.460.2230.tsv', stringsAsFactors = F)
  fimo <- filter(fimo, motif_alt_id!='')
  lapply(split(fimo,as.vector(fimo$motif_alt_id)), function(x){unique(as.vector(x$sequence_name))})
}

############################################################# map to matrix

relation.matrix <- function(a2b) {
  names(a2b) -> a
  Reduce(union,a2b) -> b
  
  lapply(a2b, function(t){
    b%in%t %>% as.numeric()
  }) -> a2b.loc
  a2b.m <- matrix(unlist(a2b.loc), ncol = length(a2b.loc))
  colnames(a2b.m) <- a
  rownames(a2b.m) <- b
  a2b.m
}
lncRNA.2.PCG <- get.lncRNA.2.PCG(0.3)
tf.2.PCG <- get.tf.2.PCG.from.fimo()
barplot(sapply(tf.2.PCG, function(x){length(x)})%>%sort)

tf.2.PCG.m <- relation.matrix(tf.2.PCG)
lncRNA.2.PCG.m <- relation.matrix(lncRNA.2.PCG)
# 
heatmap(tf.2.PCG.m)
# heatmap(lncRNA.2.PCG.m)
####################################### Fix tf.2.PCG.m & lncRNA.2.PCG.m
fix.matrix <- function(lncRNA.2.PCG.m, tf.2.PCG.m) {
  lncRNA <- colnames(lncRNA.2.PCG.m)
  lncRNA.pcg <- rownames(lncRNA.2.PCG.m)
  
  tfs <- colnames(tf.2.PCG.m)
  tf.pcg <- rownames(tf.2.PCG.m)
  
  lncRNA.2.PCG.mf <- lncRNA.2.PCG.m
  tf.2.PCG.mf <- tf.2.PCG.m[match(lncRNA.pcg, tf.pcg),]
  rownames(tf.2.PCG.mf) <- rownames(lncRNA.2.PCG.mf)
  tf.2.PCG.mf[is.na(tf.2.PCG.mf)] <- 0
  
  list(lncRNA=lncRNA.2.PCG.mf, tf=tf.2.PCG.mf)
}


dim(lncRNA.2.PCG.m)
dim(tf.2.PCG.m)
fix <- fix.matrix(lncRNA.2.PCG.m, tf.2.PCG.m)
dim(fix$lncRNA)
dim(fix$tf)
heatmap(fix$tf)
heatmap(fix$tf[,sample(ncol(fix$tf), 50)])

heatmap(fix$lncRNA)
random.lncRNA <- fix$lncRNA[,sample(ncol(fix$lncRNA), 50)]
heatmap(random.lncRNA)
# ###################################### lncRNA 2 TF Phi ###################################################
#           tf
#           1   0           
# lncRNA 1  11  10
#        0  1   0
#

phiContingency <- function(lncRNA, tf) {
  ls <- lncRNA.2.PCG.mf[,lncRNA]
  ts <- tf.2.PCG.mf[,tf]
  mt <- as.data.frame(table(ls*10+ts))
  vt <- c(
    filter(mt, Var1==11)[1,2],
    filter(mt, Var1==10)[1,2],
    filter(mt, Var1==1)[1,2],
    filter(mt, Var1==0)[1,2]
  )
  vt[is.na(vt)] <- 0
  matrix(vt, ncol = 2, byrow = T, dimnames = list(c('l1','l0'),c('t1', 't0')))
}

lncTF <- function(lncRNA.2.PCG.mf, tf.2.PCG.mf) {
  phiContingencyMulti <- function(label) {
    outer(lncRNA, tfs, function(a,b) {
      mapply(function(x,y) {
        sum((lncRNA.2.PCG.mf[,x]*10+tf.2.PCG.mf[,y])==label)
      }, a, b)
    }) -> tmp
    dimnames(tmp) <- list(lncRNA, tfs)
    tmp
  }
  
  lncRNA <- colnames(lncRNA.2.PCG.mf)
  tfs <- colnames(tf.2.PCG.mf)
  message('1. Calculate phi ...')
  outer(lncRNA, tfs, function(a,b) {
    mapply(function(x,y) {
      cor.test(lncRNA.2.PCG.mf[,x],tf.2.PCG.mf[,y])$estimate
    }, a, b)
  }) -> lncRNA.tf.phi
  dimnames(lncRNA.tf.phi) <- list(lncRNA, tfs)
  
  message('2. Calculate phi p.value ...')
  outer(lncRNA, tfs, function(a,b) {
    mapply(function(x,y) {
      cor.test(lncRNA.2.PCG.mf[,x],tf.2.PCG.mf[,y])$p.value
    }, a, b)
  }) -> lncRNA.tf.phi.p
  dimnames(lncRNA.tf.phi.p) <- list(lncRNA, tfs)
  
  message('3. Calculate  contingency table ...')
  lncRNA.tf.11 <- phiContingencyMulti(11)
  lncRNA.tf.10 <- phiContingencyMulti(10)
  lncRNA.tf.1 <- phiContingencyMulti(1)
  lncRNA.tf.0 <- phiContingencyMulti(0)
  
  #holm", "hochberg", "hommel", "bonferroni", "BH", "BY",
  #   "fdr", "none")
  lncRNA.tf.list <- reshape2::melt(lncRNA.tf.phi)
  colnames(lncRNA.tf.list) <- c('lncRNA', 'tf', 'phi')
  lncRNA.tf.list <- data.frame(lncRNA.tf.list,
                               p.value=reshape2::melt(lncRNA.tf.phi.p)$value,
                               # p.adjust=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value),
                               # p.holm=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'holm'),
                               # p.hochberg=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'hochberg'),
                               # p.hommel=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'hommel'),
                               # p.bonferroni=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'bonferroni'),
                               p.adjust=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'BH'),
                               # p.BY=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'BY'),
                               # p.fdr=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'fdr'),
                               # p.none=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'none'),
                               '11'=reshape2::melt(lncRNA.tf.11)$value,
                               '10'=reshape2::melt(lncRNA.tf.10)$value,
                               '1'=reshape2::melt(lncRNA.tf.1)$value,
                               '0'=reshape2::melt(lncRNA.tf.0)$value)
  lncRNA.tf.list
}
# lncRNA.tf.list <- lncTF(fix$lncRNA, fix$tf)
# 
# write.csv(lncRNA.tf.list, file = './data/lncRNA.tf.list.csv')


#################################################

lncTF.all <- function(s=0, tf='enricher') {
  if (tf=='enricher') {
    tf.2.PCG <- get.tf.2.PCG.from.enricher()
  } else {
    tf.2.PCG <- get.tf.2.PCG.from.fimo()
  }
  lncRNA.2.PCG <- get.lncRNA.2.PCG(s)
  lncRNA.2.PCG.m <- relation.matrix(lncRNA.2.PCG)
  tf.2.PCG.m <- relation.matrix(tf.2.PCG)
  fix <- fix.matrix(lncRNA.2.PCG.m, tf.2.PCG.m)
  print(ncol(fix$lncRNA))
  print(ncol(fix$tf))
  lncRNA.tf.list <- lncTF(fix$lncRNA, fix$tf)
  lncRNA.tf.list <- arrange(lncRNA.tf.list, p.value)
  write.csv(lncRNA.tf.list, file = paste0('./data/lncRNA.tf.list/lncRNA.tf.',tf,'.',s,'.csv'))
  lncRNA.tf.list
}

# lncTf.0.9 <- lncTF.all(0.9)
# lncTf.0.8 <- lncTF.all(0.8)
# lncTf.0.7 <- lncTF.all(0.7)
# lncTf.0.6 <- lncTF.all(0.6)
# lncTf.0.5 <- lncTF.all(0.5)
# lncTf.0.4 <- lncTF.all(0.4)
# lncTf.0.3 <- lncTF.all(0.3)
# lncTf.0.2 <- lncTF.all(0.2)
# lncTf.0.1 <- lncTF.all(0.1)
# 
# system.time(lncTf.fimo.0.9 <- lncTF.all(0.9, tf='fimo'))
# system.time(lncTf.fimo.0.8 <- lncTF.all(0.8, tf='fimo'))
# system.time(lncTf.fimo.0.7 <- lncTF.all(0.7, tf='fimo'))
# system.time(lncTf.fimo.0.6 <- lncTF.all(0.6, tf='fimo'))
system.time(lncTf.fimo.0.5 <- lncTF.all(0.5, tf='fimo'))
# system.time(lncTf.fimo.0.4 <- lncTF.all(0.4, tf='fimo'))
# system.time(lncTf.fimo.0.3 <- lncTF.all(0.3, tf='fimo'))
# system.time(lncTf.fimo.0.2 <- lncTF.all(0.2, tf='fimo'))
# system.time(lncTf.fimo.0.1 <- lncTF.all(0.1, tf='fimo'))
