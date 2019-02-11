library(pcProfile)
library(gplots)
library(pheatmap)
library(limma)
source('./R/lib.R')
#######################################################################
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
get.p <- function(x) {
  pf(x$fstatistic[1L], x$fstatistic[2L], x$fstatistic[3L], lower.tail = FALSE)
}
lmTFA <- function(gene2tf, gene2sample, cache=NULL){
  if (!is.null(cache) && file.exists(cache)) {
    readRDS(cache)
  } else {
    model <- lapply(1:ncol(gene2sample), function(i){
      x <- cbind(gene2sample[,i], gene2tf)
      colnames(x)[1] <- "EXP"
      l<-lm(EXP~.,data=as.data.frame(x))
      if (get.p(summary(l))<0.01) {
        cat(paste(i,': Model passed\n'))
        coef(summary(l))
      } else {
        cat(paste(i,': Model Failed\n'))
        NA
      }
    })
    names(model) <- colnames(gene2sample)
    if (!is.null(cache)) {
      saveRDS(model, cache)
    }
    model
  }
}

model2matrix <- function(model) {
  model.rm.na <- model[!is.na(model)]
  value <- sapply(model.rm.na, function(x){
    f <- x[,1]
    f[x[,4] > 0.05] <- 0
    f
  })
  p<-sapply(model.rm.na, function(x){
    x[,4]
  })
  value[rowSums(p<0.05)>=quantile(rowSums(p<0.05),0.95),][-1,]
}

mergeTFAM <- function(model.1, model.2) {
  m1 <- model2matrix(model.1)
  m2 <- model2matrix(model.2)
  tfs<-union(rownames(m1), rownames(m2))
  m<-cbind(m1[match(tfs, rownames(m1)),],m2[match(tfs, rownames(m2)),])
  m[is.na(m)]<-0
  rownames(m) <-tfs
  m
}

model2matrix(model.gtrd.n)->m.gtrd.n
model2matrix(model.gtrd.t)->m.gtrd.t
tf <- rownames(m.gtrd)
tf.n <- rownames(m.gtrd.n)
tf.t <- rownames(m.gtrd.t)
tf.inter <- intersect(rownames(m.gtrd.n), rownames(m.gtrd.t))
tf.t.only <- setdiff(tf, tf.n)
tf.n.only <- setdiff(tf, tf.t)

diff <- pf.get.diff()
genes <- diff$GeneID
gene2sample.m <- pf.filter.zfpkm(genes)
gene2tf.gtrd.m <- get.gene2tf.matrix(tf2gene.gtrd, genes)

model.gtrd.n <- lmTFA(gene2tf.gtrd.m, gene2sample.m[,500:551], cache = './data/tfa.model.gtrd.n.rds')
model.gtrd.t <- lmTFA(gene2tf.gtrd.m, gene2sample.m[,1:499], cache = './data/tfa.model.gtrd.t.rds')
m.gtrd <- mergeTFAM(model.gtrd.t, model.gtrd.n)
# saveRDS(m.gtrd,'./data/tfa.m.gtrd.rds')
#############################################3
m <- m.gtrd
samples <- colnames(m)

normal <- c()
tumor <- c()
for(i in samples) {
  if(grepl('11A', i)==1) {
    normal <- c(normal, i)
  } else {
    tumor <- c(tumor, i)
  }
}
normal <- normal
tumor <- tumor
nn<- length(normal)
nt<- length(tumor)
n <- nn + nt
s <- c(tumor[1:nt], normal)
annotation_col = data.frame(Tissue = factor(c(rep('Tumor',nt), rep('Normal',nn))))
rownames(annotation_col) = s
pheatmap(m[,s], 
         cluster_cols=FALSE,
         scale = 'row',
         filename= './reports/thesis/tf.activity.png',
         # width = 8,
         height = 10,
         annotation_col = annotation_col)


####################################################
################ Activity Chage ###################
##################################################
m <- m.gtrd
samples <- colnames(m)

normal <- c()
tumor <- c()
for(i in samples) {
  if(grepl('11A', i)==1) {
    normal <- c(normal, i)
  } else {
    tumor <- c(tumor, i)
  }
}
nn<- length(normal)
nt<- length(tumor)
n <- nn + nt

design <- model.matrix(~factor(c(rep('T',nt), rep('N',nn))))
fit <- lmFit(m,design)
fit <- contrasts.fit(fit,coef=2)
fit <- eBayes(fit)
ac<-topTable(fit,n=10000)
saveRDS(ac,'./data/ac.gtrd.rds')

mn <- abs(m)
design <- model.matrix(~factor(c(rep('T',nt), rep('N',nn))))
fit <- lmFit(mn,design)
fit <- contrasts.fit(fit,coef=2)
fit <- eBayes(fit)
ac2 <- topTable(fit,n=10000)
# abs
mergeAC <- function(ac) {
  a<-t(sapply(rownames(m.gtrd), function(x) {
    logFC=ac[x,'logFC']
    FDR=ac[x,'adj.P.Val']
    normal=ifelse(x%in%tf.n, '+', '-')
    tumor=ifelse(x%in%tf.t, '+', '-')
    c(x,logFC, FDR, normal, tumor)
  }))
  t <- data.frame(TF=a[,1],
                  logFC=round(as.numeric(a[,2]),3),
                  FDR=as.numeric(a[,3]),
                  Normal=a[,4],
                  Tumor=a[,5])
  t<- dplyr::arrange(t, abs(FDR))
  t$FDR <- round(t$FDR,3)
  t
}
t <- mergeAC(ac)
# write.csv(t, './reports/thesis/tfa.diff.55.csv', row.names = FALSE)
t2 <- mergeAC(ac2)
# write.csv(t2, './reports/thesis/tfa.diff.abs.55.csv', row.names = FALSE)

###################################################
diff <- pf.get.diff()

tf[pf.symbol2emsembl(tf)%in%diff$GeneID]

