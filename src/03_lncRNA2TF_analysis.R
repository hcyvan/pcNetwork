library(dplyr)
source('./lib/globals.R')
source('./lib/helpers.R')


############################################################# Code from 03_lncRNA2TF.R
load('./cache/diff.qlf.2877.pairs.rda') # diff.cor.pairs
get.lncRNA.2.PCG <- function(s=0) {
  env = globalenv()
  key = paste0('.lncRNA.2.PCG.',s)
  if (exists(key, envir = env)) {
    get(key,envir = env)
  } else {
    data <- filter(diff.cor.pairs, type1%in%config$lncRNA & type2%in%config$PCGs&significant&abs(r)>=s)
    lncRNA.2.PCG <- lapply(split(data,as.vector(data$Var1)), function(x){as.vector(x$Var2)})
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
  lncRNA.tf <- lncRNA.tf.fimo.s%>%filter(p.adjust<0.05, phi>0, X11+X10>(X11+X10+X1+X0)/10)
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
    barplot(sort(lnctp$detail$X11), main=s)
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


