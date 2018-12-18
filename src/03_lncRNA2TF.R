library(dplyr)
source('./lib/globals.R')
source('./lib/helpers.R')

############################################### GET x.2.y
get.lncRNA.2.PCG <- function(s=0, fdr=0.05) {
  env = globalenv()
  key = paste0('.lncRNA.2.PCG.',s)
  if (exists(key, envir = env)) {
    get(key,envir = env)
  } else {
    cor.pairs <- readRDS('./cache/cor.pairs.info.rds')
    data <- filter(cor.pairs, type1%in%config$lncRNA & type2%in%config$PCGs & FDR<fdr & abs(r)>=s)
    lncRNA.2.PCG <- lapply(split(data,as.vector(data$v1)), function(x){as.vector(x$v2)})
    assign(key, lncRNA.2.PCG, envir = env)
    lncRNA.2.PCG
  }
}
get.tf.2.PCG.from.fimo <- function() {
  env = globalenv()
  key = '.tf.2.PCG'
  if (exists(key, envir = env)) {
    get(key, envir = env)
  } else {
    fimo <- read.delim('./data/enricher/fimo.460.2230.tsv', stringsAsFactors = F)
    fimo <- filter(fimo, motif_alt_id!='')
    tf.2.PCG <- lapply(split(fimo,as.vector(fimo$motif_alt_id)), function(x){unique(as.vector(x$sequence_name))})
    assign(key, tf.2.PCG, envir = env)
    tf.2.PCG
  }
}

plot.lncRNA.2.PCG <- function() {
  par(mfrow=c(3,3))
  for (s in seq(0.1, 0.9, 0.1)) {
    sapply(get.lncRNA.2.PCG(s), function(x){length(x)})%>%sort->tmp
    tmp%>%barplot(main=s)
    summary(tmp)
  }
  par(mfrow=c(1,1))
}

plot.tf.2.PCG <- function() {
  sapply(get.tf.2.PCG.from.fimo(), function(x){length(x)})%>%sort->tmp
  tmp%>%barplot()
  summary(tmp)
}
# plot.lncRNA.2.PCG()
# plot.tf.2.PCG()
################################################# tf-pcg
# get.tf.2.PCG.from.enricher <- function() {
#   tf.enricher <- read.delim('./data/enricher.all.bg0.new.csv', sep = ',', stringsAsFactors = F)
#   tf.enricher.list <- setNames(split(tf.enricher, seq(nrow(tf.enricher))), tf.enricher$detail.ID)
#   
#   load('./cache/biomart.symbol.biotype.rda')
#   
#   lapply(tf.enricher.list, function(tf){
#     symbol <- str_split(tf$detail.geneID, '/')[[1]]
#     filter(biomart.symbol.biotype, hgnc_symbol%in%symbol)$ensembl_gene_id
#   })
# }

#############################################################
devtools::load_all('./package/x2y/')

lncRNA.2.PCG <- get.lncRNA.2.PCG(0.3)
tf.2.PCG <- get.tf.2.PCG.from.fimo()

tf.2.PCG.m <- x2yMatrix(tf.2.PCG)
lncRNA.2.PCG.m <- x2yMatrix(lncRNA.2.PCG)
# 
# heatmap(tf.2.PCG.m)
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
heatmap(fix$tf[,sample(ncol(fix$tf), 50)])
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

  lncRNA.tf.list <- reshape2::melt(lncRNA.tf.phi)
  colnames(lncRNA.tf.list) <- c('lncRNA', 'tf', 'phi')
  lncRNA.tf.list <- data.frame(lncRNA.tf.list,
                               p.value=reshape2::melt(lncRNA.tf.phi.p)$value,
                               FDR=p.adjust(reshape2::melt(lncRNA.tf.phi.p)$value, method = 'BH'),
                               'v11'=reshape2::melt(lncRNA.tf.11)$value,
                               'v10'=reshape2::melt(lncRNA.tf.10)$value,
                               'v1'=reshape2::melt(lncRNA.tf.1)$value,
                               'v0'=reshape2::melt(lncRNA.tf.0)$value)
  arrange(lncRNA.tf.list, FDR)
}

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
  message(paste('Calculate Phi of',ncol(fix$lncRNA), 'lncRNA and', ncol(fix$tf),'TF'))
  lncRNA.tf.list <- lncTF(fix$lncRNA, fix$tf)
  write.csv(lncRNA.tf.list, file = paste0('./data/lncRNA.tf.list/lncRNA.tf.',tf,'.',s,'.csv'))
  lncRNA.tf.list
}

system.time(lncTf.fimo.0.1 <- lncTF.all(0.1, tf='fimo'))
system.time(lncTf.fimo.0.2 <- lncTF.all(0.2, tf='fimo'))
system.time(lncTf.fimo.0.3 <- lncTF.all(0.3, tf='fimo'))
system.time(lncTf.fimo.0.4 <- lncTF.all(0.4, tf='fimo'))
system.time(lncTf.fimo.0.5 <- lncTF.all(0.5, tf='fimo'))
system.time(lncTf.fimo.0.6 <- lncTF.all(0.6, tf='fimo'))
system.time(lncTf.fimo.0.7 <- lncTF.all(0.7, tf='fimo'))
system.time(lncTf.fimo.0.8 <- lncTF.all(0.8, tf='fimo'))
system.time(lncTf.fimo.0.9 <- lncTF.all(0.9, tf='fimo'))

#############################################################################################
#############################################################################################
#############################################################################################
#############################################################################################



############################################################# Code from 03_lncRNA2TF.R

get.tf.pcgs <- function(tf) {
  get.tf.2.PCG.from.fimo()[[tf]]
}
get.lncRNA.pcgs <- function(lncRNA, l.s) {
  get.lncRNA.2.PCG(l.s)[[lncRNA]]
}
get.intersect.pcg <- function(lncRNA, tf, l.s){
  intersect(get.lncRNA.2.PCG(l.s)[[lncRNA]], get.tf.2.PCG.from.fimo()[[tf]])
}
############################################################# Read
get.lncRNA.tf.fimo <- function(s) {
  read.delim(file = paste0('./data/lncRNA.tf.list/lncRNA.tf.fimo.',s,'.csv'), sep = ',')->tmp
}
# fix.lncRNA.tf.fimo <- function(s) {
#   lncRNA.tf.fimo.s <- get.lncRNA.tf.fimo(s)
#   lncRNA.tf.fimo.s%>%select(lncRNA=Var1, tf=Var2, phi=value, p.value=p.value, p.adjust=p.adjust,X11=X11,X10=X10,X1=X1,X0=X0)->fixed
#   write.csv(fixed, file = paste0('./data/lncRNA.tf.list/lncRNA.tf.fimo.',s,'.csv'), row.names = F)
# }


lncRNA.tf.fimo.0.9 <- get.lncRNA.tf.fimo(0.9)
lncRNA.tf.fimo.0.8 <- get.lncRNA.tf.fimo(0.8)
lncRNA.tf.fimo.0.7 <- get.lncRNA.tf.fimo(0.7)
lncRNA.tf.fimo.0.6 <- get.lncRNA.tf.fimo(0.6)
lncRNA.tf.fimo.0.5 <- get.lncRNA.tf.fimo(0.5)
lncRNA.tf.fimo.0.4 <- get.lncRNA.tf.fimo(0.4)
lncRNA.tf.fimo.0.3 <- get.lncRNA.tf.fimo(0.3)
lncRNA.tf.fimo.0.2 <- get.lncRNA.tf.fimo(0.2)
lncRNA.tf.fimo.0.1 <- get.lncRNA.tf.fimo(0.1)

############################################################ parse
load('./cache/biomart.symbol.biotype.rda')
lncRNA2TF.parse<- function(lncRNA.tf.fimo.s, l.s=0.5) {
  lncRNA.tf <- lncRNA.tf.fimo.s%>%filter(FDR<0.05, phi>0, v11+v10>(v11+v10+v1+v0)/10)
  if (nrow(lncRNA.tf)==0) {
    tf=vector()
    lncRNA=vector()
    detail.inter=list()
    pcg=vector()
  } else {
    tf <- sort(unique(as.vector(lncRNA.tf$tf)))
    lncRNA <- unique(as.vector(lncRNA.tf$lncRNA))
    detail.inter <- lapply(split(lncRNA.tf, seq(nrow(lncRNA.tf))), function(x){
      get.intersect.pcg(x$lncRNA,x$tf,l.s)
    })
    names(detail.inter)<-str_c(lncRNA.tf$lncRNA, lncRNA.tf$tf, sep = '-')
    pcg <- Reduce(union, detail.inter)
  }
  index <- match(lncRNA.tf$lncRNA, biomart.symbol.biotype$ensembl_gene_id)
  symbols <- biomart.symbol.biotype[index,]$hgnc_symbol
  lncRNA.tf <- data.frame(lncRNA.tf, symbol=symbols)
  result =list(
    detail=lncRNA.tf,
    detail.inter=detail.inter,
    tf=tf,
    lncRNA=lncRNA,
    pcg=pcg
  )
  class(result) <- 'lncTP'
  result
}

print.lncTP <- function(ltp){
  cat('lncRNA', 'tf', 'pcg','lncRNA-tf', sep = '\t');cat('\n')
  cat(length(ltp$lncRNA), length(ltp$tf), length(ltp$pcg),nrow(ltp$detail), sep = '\t');cat('\n')
}


lncTP.0.9 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.9, 0.9)
lncTP.0.8 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.8, 0.8)
lncTP.0.7 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.7, 0.7)
lncTP.0.6 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.6, 0.6)
lncTP.0.5 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.5, 0.5)
lncTP.0.4 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.4, 0.4)
lncTP.0.3 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.3, 0.3)
lncTP.0.2 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.2, 0.2)
lncTP.0.1 <- lncRNA2TF.parse(lncRNA.tf.fimo.0.1, 0.1)

save(lncTP.0.1, lncTP.0.2, lncTP.0.3, lncTP.0.4, lncTP.0.5, lncTP.0.6, lncTP.0.7, lncTP.0.8, lncTP.0.9, file = './cache/lncTP.0.x.rda')
load('./cache/lncTP.0.x.rda')


lncTP.0.1
lncTP.0.2
lncTP.0.3
lncTP.0.4
lncTP.0.5
lncTP.0.6
lncTP.0.7
lncTP.0.8
lncTP.0.9

par(mfrow=c(3,3))
for(s in seq(0.1,0.9,0.1)) {
  lnctp <- get(paste0('lncTP.',s))
  cat('---------------------------',s,'----------------------------------------');cat('\n')
  print(lnctp)
  print(sort(lnctp$tf))
  print(ncol(lnctp$detail))
  if (nrow(lnctp$detail) > 0) {
    barplot(sort(lnctp$detail$v11), main=s)
  }
}


write.csv(lncTP.0.3$tf, './data/lnctp.tf.csv')
write.csv(lncTP.0.3$lncRNA, './data/lnctp.lncRNA.csv')
write.csv(lncTP.0.3$detail, './data/lnctp.detail.csv')



##############################################
## tf in pcgs
pcg.all.id <- Reduce(union, lncTP.0.3$detail.inter)
pcg.all.symbol <- biomart.symbol.biotype[match(pcg.all.id, biomart.symbol.biotype$ensembl_gene_id),'hgnc_symbol']
intersect(pcg.all.symbol, lncTP.0.3$tf)

## demo
lncTP.0.3$detail%>%filter(tf=='AR')
length(get.tf.pcgs('AR'))
length(get.lncRNA.pcgs('ENSG00000237476',0.3))



get.tf.2.PCG.from.fimo()->tf.2.PCG
get.tf.2.PCG.from.enricher()->tf.2.PCG.2

barplot(sapply(tf.2.PCG, function(x){length(x)})%>%sort())
names(tf.2.PCG)%>%sort()->fimo
barplot(sapply(tf.2.PCG.2, function(x){length(x)})%>%sort())
names(tf.2.PCG.2)%>%sort()

lncRNA.2.PCG <- get.lncRNA.2.PCG(0.5)

barplot(sapply(lncRNA.2.PCG, function(x){length(x)})%>%sort())

starbase <- read.delim('./data/lncRNA_rbp.txt',stringsAsFactors = F)
starbase$RBP%>%unique()
