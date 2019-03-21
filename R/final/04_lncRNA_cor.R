source('./R/lib.R')

lncRNA<-pf.get.diff('lncRNA')
diff<-pf.get.diff()

lncRNA.zfpkm<-pf.filter.zfpkm(lncRNA$GeneID)
zfpkm<-pf.filter.zfpkm(diff$GeneID)

cor.lncRNA.diff <- pf.multicor(t(lncRNA.zfpkm), t(zfpkm), rds = './data/cor.pearson.lncRNA.2.diff5343.zfpkm.pairs.rds')
cor.lnc.diffv1<-filter(cor.lncRNA.diff,FDR<0.01,abs(r)>=0.3)
dim(cor.lnc.diffv1)
length(unique(as.vector(cor.lnc.diffv1$v1)))
length(unique(as.vector(cor.lnc.diffv1$v2)))
table(cor.lnc.diffv1$r>0)

#######################33
cor.lnc.diffv2<-filter(cor.lnc.diffv1,v2%in%sugid)
dim(cor.lnc.diffv2)
length(unique(as.vector(cor.lnc.diffv2$v1)))
length(unique(as.vector(cor.lnc.diffv2$v2)))
intersect(as.vector(cor.lnc.diffv2$v2),sugid)

a<- with(cor.lnc.diffv2,data.frame(v1=as.vector(pf.ensembl2symbol(v1)),
                                    v2=as.vector(pf.ensembl2symbol(v2)),
                                    dist=1-abs(r),
                                    sign=ifelse(r>0,'+','-'),
                                    stringsAsFactors = FALSE))
# a<-filter(a,v2!='',v1!=v2)
b<-data.frame(node=lncRNA$symbol,type='lncRNA',stringsAsFactors = FALSE)
b$type[b$node%in%intersect(lncRNA$symbol,pf.ensembl2symbol(sugid))]<-'lncRNA-surv'
# write.csv(a, 'reports/thesis/cyto/lnc.2.surv.edge.csv',row.names = FALSE,quote=FALSE)
# write.csv(b, 'reports/thesis/cyto/lnc.node.csv',row.names = FALSE,quote = FALSE)
##################### Cor + Cis
lnc.cis<-pf.get.lnc.cis(1000)
dim(lnc.cis)
lnc.cis<-mutate(lnc.cis,key=paste0(lncRNA,gene))
cor.lnc.diffv1<-mutate(cor.lnc.diffv1,key=paste0(v1,v2))
lnc.cis.cor<-filter(cor.lnc.diffv1,key%in%lnc.cis$key)
saveRDS(lnc.cis.cor, 'support/lnc.cis.cor.rds')
dim(lnc.cis.cor)
##################### Cor + Trans
lnc.trans<-pf.get.lnc.trans()
lnc.trans<-mutate(lnc.trans,key=paste0(lncRNA,gene))
dim(lnc.trans)
length(unique(as.vector(lnc.trans$lncRNA)))
length(unique(as.vector(lnc.trans$gene)))

lnc.trans.cor<-filter(cor.lnc.diffv1,key%in%lnc.trans$key)
saveRDS(lnc.trans.cor, 'support/lnc.trans.cor.rds')
dim(lnc.trans.cor)
length(unique(as.vector(lnc.trans.cor$v1)))
length(unique(as.vector(lnc.trans.cor$v2)))
table(lnc.trans.cor$r>0)
intersect(unique(c(as.vector(lnc.trans.cor$v1),as.vector(lnc.trans.cor$v2))),sugid)

a<- with(lnc.trans.cor,data.frame(v1=as.vector(pf.ensembl2symbol(v1)),
                                   v2=as.vector(pf.ensembl2symbol(v2)),
                                   dist=1-abs(r),
                                   sign=ifelse(r>0,'+','-'),
                                   stringsAsFactors = FALSE))
a<-filter(a,v2!='',v1!=v2)
nodes<-unique(c(a$v1,a$v2))
b<-data.frame(node=nodes,
              type=ifelse(nodes%in%lncRNA$symbol,'lncRNA','gene'),
              color=ifelse(nodes%in%pf.ensembl2symbol(sugid),'surv', ifelse(nodes%in%lncRNA$symbol,'normal-lnc','normal-gene')),
              stringsAsFactor=TRUE)


write.csv(a, 'reports/thesis/cyto/lnc.2.diff.trans.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/lnc.2.diff.trans.node.csv',row.names = FALSE,quote = FALSE)
####################
cor.lnc.diffv5<-filter(cor.lnc.diffv4,v2%in%sugid|v1%in%sugid)
dim(cor.lnc.diffv5)
length(unique(as.vector(cor.lnc.diffv5$v1)))
length(unique(as.vector(cor.lnc.diffv5$v2)))

a<- with(cor.lnc.diffv5,data.frame(v1=as.vector(pf.ensembl2symbol(v1)),
                                   v2=as.vector(pf.ensembl2symbol(v2)),
                                   dist=1-abs(r),
                                   sign=ifelse(r>0,'+','-'),
                                   stringsAsFactors = FALSE))
a<-filter(a,v2!='',v1!=v2)
nodes<-unique(c(a$v1,a$v2))
b<-data.frame(node=nodes,
              type=ifelse(nodes%in%lncRNA$symbol,'lncRNA','gene'),
              color=ifelse(nodes%in%pf.ensembl2symbol(sugid),'surv', ifelse(nodes%in%lncRNA$symbol,'normal-lnc','normal-gene')),
              stringsAsFactor=TRUE)
              

write.csv(a, 'reports/thesis/cyto/lnc.2.surv.trans.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/lnc.2.surv.trans.node.csv',row.names = FALSE,quote = FALSE)

