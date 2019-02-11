library(pcProfile)
source('./lib/helpers.R')
source('./R/lib.R')


diff<- pf.get.diff()
lncRNA <- pf.get.diff('lncRNA')

data(tf2gene.gtrd)
head(tf2gene.gtrd)

cor.pairs.format <- function(s=0.3) {
    biomart <- pf.get.biomart()
    cor.pairs <- readRDS('./cache/cor.pearson.all.zfpkm.pairs.rds')
    t1 <- biomart[match(cor.pairs$v1, biomart$ensembl_gene_id), 'gene_biotype']
    t2 <- biomart[match(cor.pairs$v2, biomart$ensembl_gene_id), 'gene_biotype']
#    data.frame(cor.pairs,t1=t1,t2=t2) %>%
 #       filter((t1%in%pv.lncRNA), FDR<0.05, abs(r) >= s) %>%
  #      select(lncRNA=v1, gene=v2, r=r)
    data.frame(cor.pairs,t1=t1,t2=t2) %>%
        filter((t1%in%pv.lncRNA)) %>%
        select(lncRNA=v1, gene=v2, r=r)

}
pc <- cor.pairs.format(0)
pc.3 <- cor.pairs.format(0.3)
pc.4 <- cor.pairs.format(0.4)
pc.5 <- cor.pairs.format(0.5)
pc.6 <- cor.pairs.format(0.6)
pc.7 <- cor.pairs.format(0.7)
pc.8 <- cor.pairs.format(0.8)
pc.9 <- cor.pairs.format(0.9)


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

a3<-gene2xMatrix(x=pc.3$lncRNA,gene=pc.3$gene, value=pc.3$r, diff$GeneID)
a4<-gene2xMatrix(x=pc.4$lncRNA,gene=pc.4$gene, value=pc.4$r, diff$GeneID)
a5<-gene2xMatrix(x=pc.5$lncRNA,gene=pc.5$gene, value=pc.5$r, diff$GeneID)
a6<-gene2xMatrix(x=pc.6$lncRNA,gene=pc.6$gene, value=pc.6$r, diff$GeneID)
a7<-gene2xMatrix(x=pc.7$lncRNA,gene=pc.7$gene, value=pc.7$r, diff$GeneID)
a8<-gene2xMatrix(x=pc.8$lncRNA,gene=pc.8$gene, value=pc.8$r, diff$GeneID)
a9<-gene2xMatrix(x=pc.9$lncRNA,gene=pc.9$gene, value=pc.9$r, diff$GeneID)

a<-gene2xMatrix(x=pc$lncRNA,gene=pc$gene, value=pc$r, diff$GeneID)


b<-gene2xMatrix(x=tf2gene.gtrd$tf,gene=tf2gene.gtrd$gene, value=tf2gene.gtrd$N, diff$GeneID)


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
cs<-multicor(a,b,rds='./data/lnc2tf.cs.cor.rds',method = 'spearman')
c<-multicor(a,b,rds='./data/lnc2tf.cor.rds')

c3<-multicor(a3,b,rds='./data/lnc2tf.cor.3.rds')
c4<-multicor(a4,b,rds='./data/lnc2tf.cor.4.rds')
c5<-multicor(a5,b,rds='./data/lnc2tf.cor.5.rds')
c6<-multicor(a6,b,rds='./data/lnc2tf.cor.6.rds')
c7<-multicor(a7,b,rds='./data/lnc2tf.cor.7.rds')
c8<-multicor(a8,b,rds='./data/lnc2tf.cor.8.rds')
c9<-multicor(a9,b,rds='./data/lnc2tf.cor.9.rds')

filter(c,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
filter(cs,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))

filter(c3,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
filter(c4,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
filter(c5,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
filter(c6,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
filter(c7,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
filter(c8,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
filter(c9,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))

bb<- b[,'MYC']
aa3 <- a3[,'ENSG00000225177']
aa9 <- a9[,'ENSG00000225177']








