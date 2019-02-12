library(pcProfile)
source('./lib/helpers.R')
source('./R/lib.R')


diff<- pf.get.diff()
lncRNA <- pf.get.diff('lncRNA')

data(tf2gene.gtrd)
head(tf2gene.gtrd)

cor.pairs.format <- function(s=0.3) {
    biomart <- pf.get.biomart()
    cor.pairs <- readRDS('./data/cor.pearson.all.zfpkm.pairs.rds')
    cor.pairs <- filter(cor.pairs, FDR<=0.01, abs(r)>=s)
    t1 <- biomart[match(cor.pairs$v1, biomart$ensembl_gene_id), 'gene_biotype']
    t2 <- biomart[match(cor.pairs$v2, biomart$ensembl_gene_id), 'gene_biotype']
    tt <- data.frame(cor.pairs,t1=t1,t2=t2)
    d1 <- filter(tt, t1%in%pv.lncRNA)
    d1 <- select(d1, lncRNA=v1, gene=v2, r=r, fdr=FDR)
    d2 <- filter(tt, t2%in%pv.lncRNA, !t1%in%pv.lncRNA)
    d2 <- select(d2, lncRNA=v2, gene=v1, r=r, fdr=FDR)
    rbind(d1,d2)
}
gene2xMatrix <- function(x, gene, value, g.filter, x.filter=NULL) {
  x2gene <- data.frame(x=x,gene=gene,value=value)
  gene2x.split <- split(x2gene, x2gene$gene)[g.filter]
  xs <- Reduce(union, lapply(gene2x.split, function(x){x$x}))
  if (!is.null(x.filter)) {
    xs <- x.filter
  }
  gene2x.l <- lapply(gene2x.split, function(e){
    tmp <- e$value[match(xs,e$x)]
    if(length(tmp)==0){return(NA)}
    tmp
  })
  gene2x.m <-do.call(rbind, gene2x.l)
  rownames(gene2x.m) <- g.filter
  colnames(gene2x.m) <- xs
  gene2x.m[is.na(gene2x.m)] <- 0
  gene2x.m
}
multicor <- function(m1, m2=NULL, method= c('pearson', 'kendall', 'spearman'), rds=NA, rewrite=FALSE, verbose=TRUE) {
  if(!is.na(rds)) {
    if(file.exists(rds) && !rewrite) {
      message(paste('Load data from', rds))
      return(readRDS(rds))
    }
  }
  method <- match.arg(method)
  if (is.null(m2)) {
      data1 <- data2 <- as.matrix(m1)
  } else {
      data1 <- as.matrix(m2)
      data2 <- as.matrix(m1)
  }
  r <- matrix(nrow=ncol(data2),ncol=ncol(data1))
  if (verbose) {
      message(paste('Calculating:', nrow(data),'rounds needed!'))
      total <- ncol(data1)
  }
  i <- 0
  p <- apply(data1,2,function(x){
    i<<-i+1
    t0 <- Sys.time()
    if (verbose) {
      cat(paste0(i,'/',total,'\n'))
    }
    j<-0
    apply(data2,2,function(y){
      j<<-j+1
      if (!is.null(data2)|i > j) {
        test <- cor.test(x,y,method = method)
        r[j,i] <<- test$estimate
        test$p.value
      } else {
        NaN
      }
    })
  })
  dimnames(r) <- dimnames(p)
  list(r,p)
  r.melt <- reshape2::melt(r) %>% filter(!is.na(value))
  p.melt <- reshape2::melt(p) %>% filter(!is.na(value))
  ret <- data.frame(v1=r.melt$Var1, v2=r.melt$Var2, r=r.melt$value, p.value=p.melt$value, FDR=p.adjust(p.melt$value, method = 'BH'))
  if (!is.na(rds)) {
    saveRDS(ret, rds)
  }
  ret
}

pc3 <- cor.pairs.format(0.3)
pc4 <- cor.pairs.format(0.4)
pc5 <- cor.pairs.format(0.5)
pc6 <- cor.pairs.format(0.6)
pc7 <- cor.pairs.format(0.7)
pc8 <- cor.pairs.format(0.8)
pc9 <- cor.pairs.format(0.9)

a3<-gene2xMatrix(x=pc3$lncRNA,gene=pc3$gene, value=pc3$r, diff$GeneID)
a4<-gene2xMatrix(x=pc4$lncRNA,gene=pc4$gene, value=pc4$r, diff$GeneID)
a5<-gene2xMatrix(x=pc5$lncRNA,gene=pc5$gene, value=pc5$r, diff$GeneID)
a6<-gene2xMatrix(x=pc6$lncRNA,gene=pc6$gene, value=pc6$r, diff$GeneID)
a7<-gene2xMatrix(x=pc7$lncRNA,gene=pc7$gene, value=pc7$r, diff$GeneID)
a8<-gene2xMatrix(x=pc8$lncRNA,gene=pc8$gene, value=pc8$r, diff$GeneID)
a9<-gene2xMatrix(x=pc9$lncRNA,gene=pc9$gene, value=pc9$r, diff$GeneID)

b<-gene2xMatrix(x=tf2gene.gtrd$tf,gene=tf2gene.gtrd$gene, value=tf2gene.gtrd$N, diff$GeneID)

cs3<-multicor(a3,b,rds='./data/lnc2tf.cs.cor.0.3.rds',method = 'spearman')
cs4<-multicor(a4,b,rds='./data/lnc2tf.cs.cor.0.4.rds',method = 'spearman')
cs5<-multicor(a5,b,rds='./data/lnc2tf.cs.cor.0.5.rds',method = 'spearman')
cs6<-multicor(a6,b,rds='./data/lnc2tf.cs.cor.0.6.rds',method = 'spearman')
cs7<-multicor(a7,b,rds='./data/lnc2tf.cs.cor.0.7.rds',method = 'spearman')
cs8<-multicor(a8,b,rds='./data/lnc2tf.cs.cor.0.8.rds',method = 'spearman')
cs9<-multicor(a9,b,rds='./data/lnc2tf.cs.cor.0.9.rds',method = 'spearman')

####################### TF lncRNA cor
tfs <- unique(tf2gene.gtrd$tf)
tfs.symbol<-pf.symbol2emsembl(tfs)
tf.zfpkm<-na.omit(pf.filter.zfpkm(tfs.symbol))
rownames(tf.zfpkm) <- pf.ensembl2symbol(rownames(tf.zfpkm))
lncRNA.zfpkm <- pf.filter.zfpkm(lncRNA$GeneID)

cor.lnc2tf <- multicor(t(lncRNA.zfpkm),t(tf.zfpkm),rds='./data/cor.lnc2tf.rds')
cor.lnc2tf <- filter(cor.lnc2tf, FDR < 0.05)%>%mutate(key=paste0(v1,v2))
#######################

cs3 <- filter(cs3,FDR<=0.05)%>%mutate(key=paste0(v1,v2))
final <- filter(cs3,!key%in%cor.lnc2tf$key)
filter(final,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))


final.n<- sapply(split(final, as.vector(final$v2)),function(x){
    nrow(x)
})

cs3.n<-sapply(split(cs3, as.vector(cs3$v2)),function(x){
    nrow(x)
})
