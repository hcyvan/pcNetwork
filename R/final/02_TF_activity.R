library(pcProfile)
library(ggpubr)
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

model2matrix <- function(model,all=FALSE) {
  model.rm.na <- model[!is.na(model)]
  value <- sapply(model.rm.na, function(x){
    f <- x[,1]
    f[x[,4] > 0.05] <- 0
    f
  })
  p<-sapply(model.rm.na, function(x){
    x[,4]
  })
  if (all) {
    value[-1,]
  }else {
    value[rowSums(p<0.05)>=quantile(rowSums(p<0.05),0.75),][-1,]
  }
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
####################### TF Activty
diff<-pf.get.diff()
gene2sample.m <- pf.filter.zfpkm(diff$GeneID)
gene2tf.gtrd.m <- get.gene2tf.matrix(tf2gene.gtrd, diff$GeneID)

model.gtrd.n <- lmTFA(gene2tf.gtrd.m, gene2sample.m[,500:551], cache = './data/tfa2.model.gtrd.n.rds')
model.gtrd.t <- lmTFA(gene2tf.gtrd.m, gene2sample.m[,1:499], cache = './data/tfa2.model.gtrd.t.rds')
m.gtrd.all<-model2matrix(c(model.gtrd.n,model.gtrd.t))
tfs<-rownames(m.gtrd.all)
tfs[tfs=="`NKX2-3`"]<-"NKX2-3"
rownames(m.gtrd.all)<-tfs
######################

heatmap.2(m.gtrd.all,
          key.title = NA,
          key.xlab = 'Z-Score of signal',
          key.ylab = NA,
          trace = 'none',
          scale = 'row',
          srtCol=45,
          cexCol = 1.1,
          # adjCol = c(0.5,0.5),
          lhei = c(1,4))

###################### TF Activty Change
split.sample<- function(samples) {
  normal <- c()
  tumor <- c()
  for(i in samples) {
    if(grepl('11A', i)==1) {
      normal <- c(normal, i)
    } else {
      tumor <- c(tumor, i)
    }
  }
  list(n=normal, t=tumor)
}


m <- m.gtrd.all
samples<-split.sample(colnames(m))

nn<- length(samples$n)
nt<- length(samples$t)
n <- nn + nt

design <- model.matrix(~factor(c(rep('N',nn),rep('T',nt))))
fit <- lmFit(m,design)
fit <- contrasts.fit(fit,coef=2)
fit <- eBayes(fit)
ac<-topTable(fit,n=10000,genelist = rownames(m))
ac<- filter(ac, adj.P.Val<0.01)%>%arrange(desc(abs(logFC)))
saveRDS(ac,'data/tf.ac.rds')
ac.thesis<-dplyr::select(ac,ID,logFC)%>%mutate(logFC=round(logFC,3))
write.csv(ac.thesis,'./reports/thesis/tf.ac.thesis.csv',row.names = FALSE)
#############################################3
m <- m.gtrd.all[match(ac$ID,rownames(m.gtrd.all)),]
samples<-split.sample(colnames(m))

nn<- length(samples$n)
nt<- length(samples$t)
s<-c(samples$t,samples$n)
annotation_col = data.frame(Tissue = factor(c(rep('Tumor',nt), rep('Normal',nn))))
rownames(annotation_col) = s
pheatmap(m[,s], 
         cluster_cols=FALSE,
         scale = 'row',
         filename= './reports/thesis/tf.activity2.pdf',
         width = 8,
         height = 16,
         annotation_col = annotation_col)
##########################################3

diff.tf<-intersect(diff$GeneID,pf.symbol2emsembl(ac$ID))
pf.ensembl2symbol(diff.tf)
fpkm<-pf.filter.fpkm(pf.symbol2emsembl(ac$ID))
#--------------------
mean.fpkm<-sort(rowMeans(fpkm))
ggplot(data.frame(fpkm=mean.fpkm), aes(x=fpkm)) +
  geom_histogram(binwidth = 1) +
  labs(x='FPKM', y='TF numbers') +theme_light()+
  theme(axis.text=element_text(size=rel(1.2)),
        legend.text= element_text(size=rel(1)),
        legend.title=element_blank(),
        axis.title=element_text(size=rel(1.2)))

pcg<-pf.get.diff('pcg')
fpkm2<- pf.filter.fpkm(pcg$GeneID)
mean.fpkm2<-sort(rowMeans(fpkm2))
ggplot(data.frame(fpkm=mean.fpkm2), aes(x=fpkm)) +
  geom_histogram(binwidth = 1) +
  labs(x='FPKM', y='TF numbers') +theme_light()+
  theme(axis.text=element_text(size=rel(1.2)),
        legend.text= element_text(size=rel(1)),
        legend.title=element_blank(),
        axis.title=element_text(size=rel(1.2)))

#----------------------
logFPKM<-pf.filter.logFpkm(pf.symbol2emsembl(ac$ID))
rownames(logFPKM)<-ac$ID
get.tf<-function(tf){
  data.frame(exp=pf.filter.logFpkm(pf.symbol2emsembl(tf)),
             gene=tf,
             tissue=c(rep('Tumor',499),rep('Normal',52)),
             stringsAsFactors = FALSE)
}
tfs<-c('E2F7','HOXA4','SOX2','AR','MYC')
data<-do.call(rbind,lapply(tfs, function(x){
  get.tf(x)
}))
data<-filter(data,exp>-10)

win.metafile(filename="./reports/thesis/tf.ac.expression.emf",width=10,height=6)
ggboxplot(data,x = "gene", y = "exp",color = "tissue",
          ylab = 'log2(FPKM)', xlab = '',
          add = "jitter", add.params = list(fill = "white"),ggtheme = theme_light())+stat_compare_means(aes(group = tissue),label = "p.signif")+
theme(axis.text=element_text(size=rel(1.2)),
      legend.text= element_text(size=rel(1)),
      legend.title=element_blank(),
      axis.title=element_text(size=rel(1.2)))
dev.off()
########################################
