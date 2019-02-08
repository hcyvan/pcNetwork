
rm(list=ls());gc()
library(pcProfile)
source('./R/lib.R')

diff <- pf.get.diff()
genes <- diff$GeneID

get.gene2tf.matrix <- function(tf2gene, genes, tfs.filter=NULL){
  gene2tf.split <- split(tf2gene, tf2gene$gene)[genes]
  tfs <- Reduce(union, lapply(gene2tf.split, function(x){x$tf}))
  if (!is.null(tfs.filter)) {
    tfs <- tfs.filter
  }
  gene2tf.l <- lapply(gene2tf.split, function(x){
    tmp <- x$N[match(tfs,x$tf)]
    if(length(tmp)==0){return(NA)}
    tmp
  })
  gene2tf.m <-do.call(rbind, gene2tf.l)
  rownames(gene2tf.m) <- genes
  colnames(gene2tf.m) <- tfs
  gene2tf.m[is.na(gene2tf.m)] <- 0
  gene2tf.m
}

gene2tf.m <- get.gene2tf.matrix(tf2gene.gtrd, genes)
gene2sample.m <- pf.filter.zfpkm(genes)

expr <- do.call(rbind,lapply(1:ncol(gene2sample.m),function(i){
  data.table(sample=i,gene=1:nrow(gene2sample.m),exp=gene2sample.m[,i])
}))
tf <- data.table(gene=1:nrow(gene2tf.m),gene2tf.m)
data <- merge(expr,tf,by=c('gene'))



#############################3
library(data.table)
gene2tf.m <- as.data.table(gene2tf.m)
fun.getsample <- function(i){
  x <- cbind(gene2sample.m[,i], gene2tf.m)
  colnames(x)[1] <- "EXP"
  x
}
fun.model <- function(x){
  lm(EXP~.,data=as.data.frame(x))
}
data.l.2 <- lapply(1:10,fun.getsample)
model.l <- lapply(data.l,fun.model)
model.p <- sapply(model.l,function(x){coef(summary(x))[,4]})


##############################

fun.getsample <- function(i){
  x <- cbind(gene2sample.m[,i], gene2tf.m)
  colnames(x)[1] <- "EXP"
  x
}
fun.model <- function(x){
  lm(EXP~.,data=as.data.frame(x))
}
data.l <- lapply(1:10,fun.getsample)
model.l <- lapply(data.l,fun.model)
model.p <- sapply(model.l,function(x){coef(summary(x))[,4]})


data2 <- do.call(rbind,data.l)
model2 <- lm(EXP~.,data=as.data.frame(data2))
hist(coef(summary(model2))[,4])
