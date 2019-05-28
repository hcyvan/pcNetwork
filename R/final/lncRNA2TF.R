


################
source('./R/lib.R')
library(pcProfile)
# gene<-pf.get.diff('pcg')
gene<-pf.get.diff()
lncRNA.zfpkm<-pf.filter.zfpkm(lncRNA$GeneID[1:5])
gene.zfpkm<-pf.filter.zfpkm(gene$GeneID)
##########s
aCor<-function(lncRNA, tf, gene, exp, tfbd) {
  
}
diff.lncRNA<-pf.get.diff('lncRNA')
diff.gene<-pf.get.diff()
lncRNA<-diff.lncRNA$GeneID
# lncRNA<-c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989')
# lncRNA<-c('ENSG00000225177')
# tf<-c('MYC','AR')
# tf<-c('MYC')
gene<-diff.gene$GeneID
exp<-t(pf.get.zfpkm())
tfbd<-tf2gene.gtrd
# tfbd<-tf2gene.jaspar
r.value<-0.4
FDR.value<-0.01

lncRNA.exp<-exp[,match(lncRNA,colnames(exp))]
gene.exp<-exp[,match(gene,colnames(exp))]
tf.exp<-exp[,match(pf.symbol2emsembl(unique(as.vector(tfbd$tf))),colnames(exp))]
tf.exp<-tf.exp[,colSums(is.na(tf.exp))==0]
lncRNA2gene.pair.org<-multicor(lncRNA.exp, gene.exp, rds='./data/cor.pearson.1148lncRNA2gene5946.zfpkm.pairs.rds')
lncRNA2tf.pair.org<-multicor(lncRNA.exp, tf.exp, rds='./data/cor.pearson.1148lncRNA2tf809.zfpkm.pairs.rds')


lncRNA2gene.pair.org<-multicor(lncRNA.exp, gene.exp, rds='./data/cor.spearman.1148lncRNA2gene5946.zfpkm.pairs.rds',method = 'spearman')
lncRNA2tf.pair<-mutate(lncRNA2tf.pair.org,key=paste0(v1,pf.ensembl2symbol(v2)),tf=pf.ensembl2symbol(v2))%>%filter(FDR<0.01)

tf2gene.pair<-dplyr::select(tfbd,x=tf,gene=gene,value=N)
gene2tf.matrix<-gene2xMatrix(tf2gene.pair,gene,x.filter = c('AR'))

lncRNA2gene.pair<-filter(lncRNA2gene.pair.org, FDR<=0.01,r>=0.4)%>%dplyr::select(x=v1,gene=v2,value=r)%>%mutate(value=abs(value))
gene2lncRNA.matrix<-gene2xMatrix(lncRNA2gene.pair, gene)
acor<-multicor(gene2lncRNA.matrix, gene2tf.matrix,method = 'spearman')
acor<-mutate(acor,key=paste0(v1,v2))%>%filter(FDR<0.05)
#########################
# acor.filter<-filter(acor,FDR<0.01,!v1%in%art)%>%arrange(FDR)
acor.filter<-filter(!key%in%unique(as.vector(lncRNA2tf.pair$key)))%>%arrange(FDR)
getp(acor)
getp(acor.filter)
getp(filter(lnctf.cor.indirect,v2=='AR'))

getp<-function(acor) {
  ARAlincRNA<-pf.get.ARAlincRNA()
  P<-sum(lncRNA%in%ARAlincRNA)
  n<-nrow(acor)
  p<-sum(as.vector(acor$v1)%in%ARAlincRNA)
  bg<-c(rep(0,1148-P),rep(1,P))
  v<-sapply(1:200000, function(x){sum(sample(bg,n))})
  Pvalue<-approxfun(density(v))(p)
  cat('ARAlincRNA: ',length(ARAlincRNA),'; ARAlincRNA in diff: ',P,'; AR2TF: ',n,'; intersect: ', p, '; P-value: ', Pvalue,'\n')
}
###################3
ARAlincRNA<-pf.get.ARAlincRNA()
ARAlincRNA.df<-data.frame(
  id=ARAlincRNA,
  symbol=pf.ensembl2symbol(ARAlincRNA),
  ActivityCorralation=ifelse(as.vector(ARAlincRNA%in%lnctf.cor.indirect$v1),1,0),
  DiffLncRNA=ifelse(as.vector(ARAlincRNA%in%lncRNA),1,0))
write.csv(ARAlincRNA.df,'./reports/thesis/ARAlincRNA.csv',row.names = FALSE)

#----
acor.filter<-filter(acor,FDR<0.01)%>%arrange(FDR)
filter(acor,v2=='AR')
filter(acor.filter,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
cor.test(gene2lncRNA.matrix[,'ENSG00000225177'],gene2tf.matrix[,'KDM5A'],method = 'spearman')


plot(gene2lncRNA.matrix[,'ENSG00000225177'],gene2tf.matrix[,'MYC'])
plot(rank(gene2lncRNA.matrix[,'ENSG00000225177']),rank(gene2tf.matrix[,'MYC']))
plot(rank(gene2lncRNA.matrix[,'ENSG00000225177']),rank(gene2tf.matrix[,'KDM5A']))
plot((gene2lncRNA.matrix[,'ENSG00000225177']),(gene2tf.matrix[,'MYC']))
plot((gene2lncRNA.matrix[,'ENSG00000270933']),(gene2tf.matrix[,'MYC']))

###########3
cor.test(gene2lncRNA.matrix[,'ENSG00000225177'],gene2tf.matrix[,'KDM5A'],method = 'spearman')
library("ggpubr")
my_data=data.frame(lncRNA=gene2lncRNA.matrix[,'ENSG00000225177'],tf=gene2tf.matrix[,'KDM5A'])
ggscatter(my_data, x = "lncRNA", y = "tf", 
          add = "reg.line", conf.int = TRUE, 
          cor.coef = TRUE, cor.method = "spearman",
          xlab = "Miles/(US) gallon", ylab = "Weight (1000 lbs)")
my_data=data.frame(lncRNA=rank(gene2lncRNA.matrix[,'ENSG00000225177']),tf=rank(gene2tf.matrix[,'MYC']))
plot()
ggplot(my_data,aes(x=lncRNA,y=tf))+geom_point()+theme_bw()
ggplot(my_data,aes(x=lncRNA,y=tf))+geom_hex(bins=50) + theme_bw()
ggplot(my_data,aes(x=lncRNA,y=tf))+geom_density_2d()




################
gene2xMatrix <- function(x2gene.pair, g.filter, x.filter=NULL) {
  gene2x.split <- split(x2gene.pair, x2gene.pair$gene)[g.filter]
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
    total <- ncol(data1)
    message(paste('Calculating:', total,'rounds needed!'))
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
      test <- cor.test(x,y,method = method)
      r[j,i] <<- test$estimate
      test$p.value
    })
  })
  dimnames(r) <- dimnames(p)
  r.melt <- reshape2::melt(r) %>% filter(!is.na(value))
  p.melt <- reshape2::melt(p) %>% filter(!is.na(value))
  ret <- data.frame(v1=r.melt$Var1, v2=r.melt$Var2, r=r.melt$value, p.value=p.melt$value, FDR=p.adjust(p.melt$value, method = 'BH'),stringsAsFactors = FALSE)
  if (!is.na(rds)) {
    saveRDS(ret, rds)
  }
  ret
}
multicor(gene2lncRNA.matrix, gene2tf.matrix,method = 'spearman')
