library(pcProfile)
library(ggplot2)
library(cowplot)
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
    # cor.pairs <- filter(cor.pairs, FDR<=0.01)
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

interNum<- function(lnc, tf, a, b) {
  ag<-a[,lnc]
  bg<-b[,tf]
  ag[ag!=0]<-1
  bg[bg!=0]<-1
  tag <-ag*10+bg
  c(sum(tag==11),sum(tag==10),sum(tag==1),sum(tag==0))
}

pc3 <- cor.pairs.format(0.3)
pc4 <- cor.pairs.format(0.4)

a3<-gene2xMatrix(x=pc3$lncRNA,gene=pc3$gene, value=pc3$r, diff$GeneID)
a4<-gene2xMatrix(x=pc4$lncRNA,gene=pc4$gene, value=pc4$r, diff$GeneID)

b<-gene2xMatrix(x=tf2gene.gtrd$tf,gene=tf2gene.gtrd$gene, value=tf2gene.gtrd$N, diff$GeneID)

cs3.abs<-multicor(abs(a3),b,rds='./data/lnc2tf.cs.cor.3.abs.rds',method = 'spearman')
cs4.abs<-multicor(abs(a4),b,rds='./data/lnc2tf.cs.cor.4.abs.rds',method = 'spearman')
cs5.abs<-multicor(abs(a5),b,rds='./data/lnc2tf.cs.cor.5.abs.rds',method = 'spearman')
cs6.abs<-multicor(abs(a6),b,rds='./data/lnc2tf.cs.cor.6.abs.rds',method = 'spearman')
cs7.abs<-multicor(abs(a7),b,rds='./data/lnc2tf.cs.cor.7.abs.rds',method = 'spearman')
cs8.abs<-multicor(abs(a8),b,rds='./data/lnc2tf.cs.cor.8.abs.rds',method = 'spearman')
cs9.abs<-multicor(abs(a9),b,rds='./data/lnc2tf.cs.cor.9.abs.rds',method = 'spearman')
####################### TF lncRNA cor
tfs <- unique(tf2gene.gtrd$tf)
tfs.symbol<-pf.symbol2emsembl(tfs)
tf.zfpkm<-na.omit(pf.filter.zfpkm(tfs.symbol))
rownames(tf.zfpkm) <- pf.ensembl2symbol(rownames(tf.zfpkm))
lncRNA.zfpkm <- pf.filter.zfpkm(lncRNA$GeneID)

cor.lnc2tf <- multicor(t(lncRNA.zfpkm),t(tf.zfpkm),rds='./data/cor.lnc2tf.rds')
cor.lnc2tf <- filter(cor.lnc2tf, FDR < 0.05)%>%mutate(key=paste0(v1,v2))
#######################
final<-mutate(cs4.abs,key=paste0(v1,v2))%>%filter(FDR<=1e-5,!key%in%cor.lnc2tf$key)
dim(final)
filter(final,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
tfs0<- sort(sapply(split(final, as.vector(final$v2)),function(x){nrow(x)}))
lncs0<- sort(sapply(split(final, as.vector(final$v1)),function(x){nrow(x)}))
length(tfs0)
length(lncs0)
tfs0['MYC']
tfs0['AR']
tfs0['TP53']



##########################################
actf<-read.delim('./reports/thesis/tfa.diff.55.csv', stringsAsFactors=FALSE,sep=',')$TF
tfs0[tfs0>quantile(tfs0,0.90)]
final9 <- filter(final,v2%in%names(tfs0[tfs0>quantile(tfs0,0.90)]))
lncs1<- sort(sapply(split(final9, as.vector(final9$v1)),function(x){nrow(x)}))
lncs1['ENSG00000225177']
lncs1['ENSG00000270933']
final95 <- filter(final9,v1%in%names(lncs1[lncs1<quantile(lncs1,0.5)]))

tfs<- sort(sapply(split(final95, as.vector(final95$v2)),function(x){nrow(x)}))
lncs<- sort(sapply(split(final95, as.vector(final95$v1)),function(x){nrow(x)}))
length(tfs)
length(lncs)
tfs['MYC']
tfs['AR']
tfs['TP53']
lncs['ENSG00000225177']
lncs['ENSG00000270933']

intersect(unique(as.vector(final95$v2)),actf)

write.csv(final95,'./reports/thesis/lnctf.cyto.csv', row.names = FALSE)

###########################33
p1<-ggplot(data.frame(tf=names(tfs0),N=tfs0), aes(x=N)) +
  geom_histogram(binwidth = 1) +
  labs(x='lncRNA Number', y='TF Number')
p2<-ggplot(data.frame(lnc=names(lncs0),N=lncs0), aes(x=N)) +
  geom_histogram(binwidth = 1) +
  labs(x='TF Number', y='lncRNA Number')

p3<-ggplot(data.frame(tf=names(tfs),N=tfs), aes(x=N)) +
  geom_histogram(binwidth = 1) +
  labs(x='lncRNA Number', y='TF Number')
p4<-ggplot(data.frame(lnc=names(lncs),N=lncs), aes(x=N)) +
  geom_histogram(binwidth = 1) +
  labs(x='TF Number', y='lncRNA Number')
win.metafile(filename="./reports/thesis//lnc2tf.dist.emf",width=7,height=7)
plot_grid(p1,p2,p3,p4,labels = c('A','B','C','D'))
dev.off()



#################3


data <- data.frame(
  Protein=c('c-Myc','IgG'),
  lncRNA=c('AC010331.1-201','AC010331.1-201','NR_033896','NR_033896','AC010719.1','AC010719.1'),
  mean=c(0.012,0.025,0.077,0.015,0.003,0.011),
  sd=c(0.00607,0.00764,0.00138,0.00312,0.0009,0.00098)
)


p<-ggplot(data, aes(x=lncRNA, y=mean, fill=Protein)) + 
  geom_bar(stat="identity", color="black", 
           position=position_dodge()) +
  geom_errorbar(aes(ymin=mean-sd, ymax=mean+sd), width=.2,
                position=position_dodge(.9)) +
  labs(y = "%Input",x='')+
  theme(axis.title.y = element_text(size = 25), axis.title=element_text(size = 20)) +
  scale_fill_manual(values=c('#E69F00','#999999'))
ggsave('./reports/thesis/RIP-anti-c-Myc-RT-PCR.emf',p)
