source('./R/lib.R')

du <- pf.get.du()
ou <- pf.get.ou()
dug <- pf.get.dug()
oug <- pf.get.oug()
diff <- pf.get.diff()
############
primers<-read.csv('data/190422/data.csv',header = FALSE)
primers$V2%in%dug$gene
primers$V2%in%oug$gene


pdu<-du[match(primers$V2,du$gene),]
pou<-ou[match(primers$V2,ou$gene),]
diff.lnc<-diff[match(primers$V2,diff$GeneID),]


candidate<-data.frame(primers, logFC=diff.lnc$logFC,Pvalue=diff.lnc$PValue,pdu$ph,pdu$ph.p, pou$ph,pou$ph.p)

write.csv(candidate,'data/190422/candidate.csv',row.names = FALSE)


for(gene in candidate$V2) {
  pf.plot.survival(gene,type='os',dir='data/190422/surv/os/',file.type = 'png')
  pf.plot.survival(gene,type='dfs',dir='data/190422/surv/dfs/',file.type = 'png')
  pf.plot.diff(gene,dir='data/190422/surv/diff/', file.type = 'png')
}

### 20190429
#### 1
wd<-read.csv('data/190422/wd.csv',header = TRUE)
gene.list <- pf.symbol2emsembl(wd$symbol)

pdu<-du[match(gene.list,du$gene),]
pou<-ou[match(gene.list,ou$gene),]
diff.lnc<-diff[match(gene.list,diff$GeneID),]

wd <- data.frame(wd, pdu$ph,pdu$ph.p, pou$ph,pou$ph.p)
write.csv(wd,'data/190422/wd_surv/wd_surv.csv',row.names = FALSE)
for(gene in gene.list) {
  pf.plot.survival(gene,type='os',dir='data/190422/wd_surv/os/',file.type = 'png')
  pf.plot.survival(gene,type='dfs',dir='data/190422/wd_surv/dfs/',file.type = 'png')
  pf.plot.diff(gene,dir='data/190422/wd_surv/diff/', file.type = 'png')
}
#### 2
read <- data("tf2gene.gtrd")
dim(filter(tf2gene.gtrd, tf=="AR"))
ar75<-quantile(filter(tf2gene.gtrd,tf=="AR")$N,0.75)
dim(filter(tf2gene.gtrd, tf=="AR", N>=ar75))

head(filter(tf2gene.gtrd, tf=="CREB3L1"))

dim(filter(tf2gene.gtrd, tf=="BACH1"))

i<-1
for(tf_name in unique(tf2gene.gtrd$tf)){
  print(i)
  tf2gene <- filter(tf2gene.gtrd, tf==tf_name)
  res <- data.frame(id=tf2gene$gene,symbol=pf.ensembl2symbol(tf2gene$gene),N=tf2gene$N)
  write.csv(res,file = file.path('data/190422/tf_chip_seq/',paste0(tf_name,'.csv')),row.names = FALSE)
  i<-i+1
}
