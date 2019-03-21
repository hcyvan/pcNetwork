library(pcProfile)
library(ggplot2)
library(cowplot)
source('./lib/helpers.R')
source('./R/lib.R')


diff<- pf.get.diff()
lncRNA <- pf.get.diff('lncRNA')
actf <- pf.get.ac()
sugid<-pf.get.sugid()

data(tf2gene.gtrd)
head(tf2gene.gtrd)
tf2gene<-filter(tf2gene.gtrd,tf%in%actf$ID)

cor.pairs.format <- function(s=0.3) {
    biomart <- pf.get.biomart()
    # cor.pairs <- readRDS('./data/cor.pearson.all.zfpkm.pairs.rds')
    cor.pairs<-readRDS('./data/cor.pearson.lncRNA.pcg.5343.zfpkm.pairs.rds')
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
pc5 <- cor.pairs.format(0.5)
pc6 <- cor.pairs.format(0.6)
pc7 <- cor.pairs.format(0.7)
pc8 <- cor.pairs.format(0.8)

a3<-gene2xMatrix(x=pc3$lncRNA,gene=pc3$gene, value=pc3$r, diff$GeneID)
a4<-gene2xMatrix(x=pc4$lncRNA,gene=pc4$gene, value=pc4$r, diff$GeneID)
a5<-gene2xMatrix(x=pc5$lncRNA,gene=pc5$gene, value=pc5$r, diff$GeneID)
a6<-gene2xMatrix(x=pc6$lncRNA,gene=pc6$gene, value=pc6$r, diff$GeneID)
a7<-gene2xMatrix(x=pc7$lncRNA,gene=pc7$gene, value=pc7$r, diff$GeneID)
a8<-gene2xMatrix(x=pc8$lncRNA,gene=pc8$gene, value=pc8$r, diff$GeneID)

b<-gene2xMatrix(x=tf2gene$tf,gene=tf2gene$gene, value=tf2gene$N, diff$GeneID)

cs3.abs<-pf.multicor(abs(a3),b,rds='./data/lnc2tf.cs5946.cor.3.abs.rds',method = 'spearman')
cs4.abs<-pf.multicor(abs(a4),b,rds='./data/lnc2actf.cs5946.cor.4.abs.rds',method = 'spearman')
cs5.abs<-pf.multicor(abs(a5),b,rds='./data/lnc2tf.cs5946.cor.5.abs.rds',method = 'spearman')
cs6.abs<-pf.multicor(abs(a6),b,rds='./data/lnc2actf.cs5946.cor.6.abs.rds',method = 'spearman')
cs7.abs<-pf.multicor(abs(a7),b,rds='./data/lnc2tf.cs5946.cor.7.abs.rds',method = 'spearman')
cs8.abs<-pf.multicor(abs(a8),b,rds='./data/lnc2tf.cs5946.cor.8.abs.rds',method = 'spearman')
cs9.abs<-multicor(abs(a9),b,rds='./data/lnc2tf.cs.cor.9.abs.rds',method = 'spearman')


filter(cs4.abs,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))

########################33
lnctfv1<-filter(cs4.abs,FDR<0.01)
dim(lnctfv1)
filter(tf2gene, tf=='MYC',gene=='ENSG00000099840')
################# lncRNA vs. TF  (Direct correlation)
actf.zfpkm<-na.omit(pf.filter.zfpkm(pf.symbol2emsembl(actf$ID)))
rownames(actf.zfpkm) <- pf.ensembl2symbol(rownames(actf.zfpkm))
lncRNA.zfpkm <- pf.filter.zfpkm(lncRNA$GeneID)
cor.lnc2actf <- pf.multicor(t(lncRNA.zfpkm),t(actf.zfpkm),rds='./data/cor.lnc2actf.rds')
cor.lnc2actf <- filter(cor.lnc2actf, FDR < 0.05)%>%mutate(key=paste0(v1,v2))
################ lncRNA vs. TF  (Indirect correlation)
lnctf.cor.indirect<-mutate(lnctfv1,key=paste0(v1,v2))%>%filter(FDR<=0.05,!key%in%cor.lnc2actf$key)
# saveRDS(lnctf.cor.indirect,'support/lnctf.cor.indirect.rds')
dim(lnctf.cor.indirect)
filter(lnctf.cor.indirect,v2=='MYC',v1%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
tfs0<- sort(sapply(split(lnctf.cor.indirect, as.vector(lnctf.cor.indirect$v2)),function(x){nrow(x)}))
lncs0<- sort(sapply(split(lnctf.cor.indirect, as.vector(lnctf.cor.indirect$v1)),function(x){nrow(x)}))
length(tfs0)
length(lncs0)
lncs0['ENSG00000225177']
tfs0['MYC']
tfs0['AR']
tfs0['TP53']
############################
lnc845.zfpkm<-na.omit(pf.filter.zfpkm(names(lncs0)))
cor.lnc2lnc <- pf.multicor(t(lnc845.zfpkm),rds='./data/cor.lnc2lnc.rds')
dim(cor.lnc2lnc)
cor.lnc2lnc <- filter(cor.lnc2lnc, FDR < 0.05)
dim(cor.lnc2lnc)
mm<-cor(t(lnc845.zfpkm))
s<-sample(1:845,100)
a<-c("ENSG00000225177","ENSG00000277383","ENSG00000270933")
filter(cor.lnc2lnc,v1%in%a,v2%in%a)

########################### Indirect + Cis
lnc.cis.cor<-pf.get.lnc.cis.cor()
head(lnctf.cor.indirect);dim(lnctf.cor.indirect)
head(lnc.cis.cor);dim(lnc.cis.cor)
head(tf2gene);dim(tf2gene)

lnc.tf.gene.cis<-left_join(lnc.cis.cor,lnctf.cor.indirect,by=c('v1'='v1'))%>%
  select(lnc=v1,tf=v2.y,gene=v2.x)%>%
  mutate(key=paste0(tf,gene))%>%
  filter(key%in%paste0(tf2gene$tf,tf2gene$gene))%>%
  select(lnc,tf,gene)
saveRDS(lnc.tf.gene.cis,'support/lnc.tf.gene.cis.rds')
# #ENSG00000227392   TAF1 ENSG00000105707
# filter(lnc.cis.cor,v1=='ENSG00000227392',v2=='ENSG00000105707')
# filter(tf2gene,tf=='TAF1',gene=='ENSG00000105707')
# filter(lnctf.cor.indirect,v1=='ENSG00000227392',v2=='TAF1')

dim(a<-distinct(lnc.tf.gene.cis,lnc,tf))
length(unique(as.vector(a$lnc)))
length(unique(as.vector(a$tf)))
dim(lnc.tf.gene.cis)

#-----------
tmp<- lnc.tf.gene.cis
tmp$lnc<-pf.ensembl2symbol(tmp$lnc)
tmp$gene<-pf.ensembl2symbol(tmp$gene)
a<-data.frame(v1=c(as.vector(tmp$lnc),
                   as.vector(tmp$tf)),
              v2=c(as.vector(tmp$tf),
                   as.vector(tmp$gene)),
              key=c(rep('b',nrow(tmp)),rep('i',nrow(tmp))),
              stringsAsFactors = FALSE)
nodes<-unique(c(a$v1,a$v2))
b<-data.frame(node=nodes,
              type=ifelse(nodes%in%lncRNA$symbol,'lncRNA',
                          ifelse(nodes%in%ac$ID,'tf','gene')),
              color=ifelse(nodes%in%pf.ensembl2symbol(sugid),'surv',
                           ifelse(nodes%in%lncRNA$symbol,'normal-lnc',
                                  ifelse(nodes%in%ac$ID,'normal-tf','normal-gene'))),
              stringsAsFactors=FALSE)
write.csv(a, 'reports/thesis/cyto/lnc.tf.gene.cis.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/lnc.tf.gene.cis.node.csv',row.names = FALSE,quote = FALSE)

#---------
lnc.tf.gene.cis.s<-mutate(lnc.tf.gene.cis,lnc=pf.ensembl2symbol(lnc),gene=pf.ensembl2symbol(gene))
sugid.s<-pf.ensembl2symbol(sugid)
lnc.tf.gene.cis.sugid<-filter(lnc.tf.gene.cis.s, lnc%in%sugid.s|tf%in%sugid.s|gene%in%sugid.s)
lnc.tf.gene.cis.sugid$tf<-as.vector(lnc.tf.gene.cis.sugid$tf)
lnc.tf.gene.cis.sugid<-mutate(lnc.tf.gene.cis.sugid,lnc=ifelse(lnc%in%sugid.s,paste0(lnc,'*'),lnc),
                          tf=ifelse(tf%in%sugid.s,paste0(tf,'*'),tf),
                          gene=ifelse(gene%in%sugid.s,paste0(gene,'*'),gene))

dim(lnc.tf.gene.cis.sugid)
length(unique(as.vector(lnc.tf.gene.cis.sugid$lnc)))
length(unique(as.vector(lnc.tf.gene.cis.sugid$tf)))
length(unique(as.vector(lnc.tf.gene.cis.sugid$gene)))
lnc.tf.gene.cis.sugid<-arrange(lnc.tf.gene.cis.sugid,tf,gene,lnc)
write.csv(lnc.tf.gene.cis.sugid,file = './reports/thesis/lnc.tf.gene.cis.sugid.csv',row.names = FALSE)

############
a<-table(filter(tf2gene.gtrd,tf=='ESR1')$N)
########################### lncRNA Trans
lnc.trans.cor<-pf.get.lnc.trans.cor()
head(lnctf.cor.indirect);dim(lnctf.cor.indirect)
head(lnc.trans.cor);dim(lnc.trans.cor)
head(tf2gene);dim(tf2gene)


lnc.tf.gene.trans<-inner_join(lnc.trans.cor,lnctf.cor.indirect,by=c('v1'='v1'))%>%
  select(lnc=v1,tf=v2.y,gene=v2.x)%>%
  mutate(key=paste0(tf,gene))%>%
  filter(key%in%paste0(tf2gene$tf,tf2gene$gene))%>%
  select(lnc,tf,gene)
saveRDS(lnc.tf.gene.trans,'support/lnc.tf.gene.trans.rds')

dim(a<-distinct(lnc.tf.gene.trans,lnc,tf))
length(unique(as.vector(a$lnc)))
length(unique(as.vector(a$tf)))
dim(lnc.tf.gene)

filter(lnc.tf.gene.trans, tf=='MYC',lnc%in%c('ENSG00000225177','ENSG00000277383','ENSG00000270933','ENSG00000197989'))
pf.ensembl2symbol('ENSG00000277383')
pf.ensembl2symbol('ENSG00000270933')

lnc.tf.gene.s<-mutate(lnc.tf.gene.trans,lnc=pf.ensembl2symbol(lnc),gene=pf.ensembl2symbol(gene))
lnc.tf.gene.s.1<-lnc.tf.gene.s
a<-data.frame(v1=c(as.vector(lnc.tf.gene.s.1$lnc),
                   as.vector(lnc.tf.gene.s.1$tf)),
              v2=c(as.vector(lnc.tf.gene.s.1$tf),
                   as.vector(lnc.tf.gene.s.1$gene)),
              key=c(rep('b',nrow(lnc.tf.gene.s.1)),rep('i',nrow(lnc.tf.gene.s.1))),
              stringsAsFactors = FALSE)
v1<-sort(table(a$v1))
v2<-sort(table(a$v2))
nodes<-unique(c(a$v1,a$v2))
b<-data.frame(node=nodes,
              type=ifelse(nodes%in%lncRNA$symbol,'lncRNA',
                          ifelse(nodes%in%ac$ID,'tf','gene')),
              color=ifelse(nodes%in%pf.ensembl2symbol(sugid),'surv',
                           ifelse(nodes%in%lncRNA$symbol,'normal-lnc',
                                  ifelse(nodes%in%ac$ID,'normal-tf','normal-gene'))),
              stringsAsFactors=FALSE)
write.csv(a, 'reports/thesis/cyto/lnc.tf.gene.trans.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/lnc.tf.gene.trans.node.csv',row.names = FALSE,quote = FALSE)
##################
filter(lnc.tf.gene.s,tf=='MYC',lnc==pf.ensembl2symbol('ENSG00000225177'))
sugid.s<-pf.ensembl2symbol(sugid)
lnc.tf.gene.sugid<-filter(lnc.tf.gene.s, lnc%in%sugid.s|tf%in%sugid.s|gene%in%sugid.s)
lnc.tf.gene.sugid$tf<-as.vector(lnc.tf.gene.sugid$tf)
lnc.tf.gene.sugid<-mutate(lnc.tf.gene.sugid,lnc=ifelse(lnc%in%sugid.s,paste0(lnc,'*'),lnc),
                          tf=ifelse(tf%in%sugid.s,paste0(tf,'*'),tf),
                          gene=ifelse(gene%in%sugid.s,paste0(gene,'*'),gene))

dim(lnc.tf.gene.sugid)
length(unique(as.vector(lnc.tf.gene.sugid$lnc)))
length(unique(as.vector(lnc.tf.gene.sugid$tf)))
length(unique(as.vector(lnc.tf.gene.sugid$gene)))
lnc.tf.gene.sugid<-arrange(lnc.tf.gene.sugid,tf,gene,lnc)
write.csv(lnc.tf.gene.sugid,file = './reports/thesis/lnc.tf.gene.sugid.csv',row.names = FALSE)



tfs<-sort(table(as.vector(lnc.tf.gene.sugid$tf)))
tfs
tfs.10<-names(tfs)[tfs>=15]
lnc.tf.gene.sugid.f<-filter(lnc.tf.gene.sugid,tf%in%tfs.10)
a<-data.frame(v1=c(as.vector(lnc.tf.gene.sugid.f$lnc),
                   as.vector(lnc.tf.gene.sugid.f$tf)),
              v2=c(as.vector(lnc.tf.gene.sugid.f$tf),
                   as.vector(lnc.tf.gene.sugid.f$gene)),
              key=c(rep('b',nrow(lnc.tf.gene.sugid.f)),rep('i',nrow(lnc.tf.gene.sugid.f))),
              stringsAsFactors = FALSE)
nodes<-unique(c(a$v1,a$v2))
b<-data.frame(node=nodes,
              type=ifelse(nodes%in%lncRNA$symbol,'lncRNA',
                          ifelse(nodes%in%ac$ID,'tf','gene')),
              color=ifelse(nodes%in%pf.ensembl2symbol(sugid),'surv',
                           ifelse(nodes%in%lncRNA$symbol,'normal-lnc',
                                  ifelse(nodes%in%ac$ID,'normal-tf','normal-gene'))),
              stringsAsFactors=FALSE)
write.csv(a, 'reports/thesis/cyto/lnc.tf.gene.surv.trans.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/lnc.tf.gene.surv.trans.node.csv',row.names = FALSE,quote = FALSE)

###############################
##############################
##########################
lnctfv3 <- filter(lnctfv2,v1%in%names(lncs0[lncs0<quantile(lncs0,0.85)]))



##########################################
# actf<-read.delim('./reports/thesis/tfa.diff.55.csv', stringsAsFactors=FALSE,sep=',')$TF
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
  lncRNA=c('AC010331.1-201','AC010331.1-201','AL590617.2','AL590617.2','AC010719.1','AC010719.1'),
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
