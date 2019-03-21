source('./R/lib.R')
library(pcProfile)

ac<- pf.get.ac()
diff<-pf.get.diff()
sugid <-pf.get.sugid()

data(tf2gene.gtrd)
head(tf2gene.gtrd)
tf2gene<-filter(tf2gene.gtrd,tf%in%ac$ID)
tf2gene<-mutate(tf2gene,key=paste0(tf,gene))

actf.zfpkm<-pf.filter.zfpkm(pf.symbol2emsembl(ac$ID))
zfpkm<-pf.filter.zfpkm(diff$GeneID)


cor.actf.diff <- pf.multicor(t(actf.zfpkm), t(zfpkm), rds = './data/cor.pearson.actf.2.diff5343.zfpkm.pairs.rds')
cor.actf.diff<-mutate(cor.actf.diff,v1=pf.ensembl2symbol(v1),key=paste0(v1,v2))
##############################################################
# cor.actf.diffv1<-filter(cor.actf.diff,FDR<0.01,abs(r)>=0.6,key%in%tf2gene$key)
cor.actf.diffv1<-filter(cor.actf.diff,FDR<0.01,abs(r)>=0.6)
dim(cor.actf.diffv1)
length(unique(as.vector(cor.actf.diffv1$v1)))
length(unique(as.vector(cor.actf.diffv1$v2)))
table(cor.actf.diffv1$r>0)

############################################################
cor.actf.diffv2<-filter(cor.actf.diffv1,v2%in%sugid|pf.symbol2emsembl(v1)%in%sugid)
dim(cor.actf.diffv2)
length(unique(as.vector(cor.actf.diffv2$v1)))
length(unique(as.vector(cor.actf.diffv2$v2)))
filter(cor.actf.diffv2,v1=='NSD2')

a<- with(cor.actf.diffv2,data.frame(v1=v1,
                                   v2=as.vector(pf.ensembl2symbol(v2)),
                                   dist=1-abs(r),
                                   sign=ifelse(r>0,'+','-'),
                                   stringsAsFactors = FALSE))
a<-filter(a,v2!='',v1!=v2)
nodes<-unique(c(a$v1,a$v2))
b<-data.frame(node=nodes,
              type=ifelse(nodes%in%tfs,'tfs','gene'),
              color=ifelse(nodes%in%pf.ensembl2symbol(sugid),'surv', ifelse(nodes%in%tfs,'normal-tfs','normal-gene')),
              stringsAsFactors=TRUE)
write.csv(a, 'reports/thesis/cyto/actf.2.surv.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/actf.2.surv.node.csv',row.names = FALSE,quote = FALSE)

################################################################
cor.actf.diffv3<-filter(cor.actf.diffv1,key%in%tf2gene$key)
dim(cor.actf.diffv3)
length(unique(as.vector(cor.actf.diffv3$v1)))
length(unique(as.vector(cor.actf.diffv3$v2)))
table(cor.actf.diffv3$r>0)

a<- with(cor.actf.diffv3,data.frame(v1=v1,
                                    v2=as.vector(pf.ensembl2symbol(v2)),
                                    dist=1-abs(r),
                                    sign=ifelse(r>0,'+','-'),
                                    stringsAsFactors = FALSE))
a<-filter(a,v2!='',v1!=v2)
nodes<-unique(c(a$v1,a$v2))
b<-data.frame(node=nodes,
              type=ifelse(nodes%in%tfs,'tfs','gene'),
              color=ifelse(nodes%in%pf.ensembl2symbol(sugid),'surv', ifelse(nodes%in%tfs,'normal-tfs','normal-gene')),
              stringsAsFactors=TRUE)
write.csv(a, 'reports/thesis/cyto/actf.direct.2.diff.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/actf.direct.2.diff.node.csv',row.names = FALSE,quote = FALSE)
##################################
cor.actf.diffv4<-filter(cor.actf.diffv3,v2%in%sugid|pf.symbol2emsembl(v1)%in%sugid)
dim(cor.actf.diffv4)
length(unique(as.vector(cor.actf.diffv4$v1)))
length(unique(as.vector(cor.actf.diffv4$v2)))
intersect(pf.ensembl2symbol(sugid),ac$ID)

a<- with(cor.actf.diffv4,data.frame(v1=v1,
                                    v2=as.vector(pf.ensembl2symbol(v2)),
                                    dist=1-abs(r),
                                    sign=ifelse(r>0,'+','-'),
                                    stringsAsFactors = FALSE))
a<-filter(a,v2!='',v1!=v2)
nodes<-unique(c(a$v1,a$v2))
b<-data.frame(node=nodes,
              type=ifelse(nodes%in%tfs,'tfs','gene'),
              color=ifelse(nodes%in%pf.ensembl2symbol(sugid),'surv', ifelse(nodes%in%tfs,'normal-tfs','normal-gene')),
              stringsAsFactors=TRUE)
write.csv(a, 'reports/thesis/cyto/actf.direct.2.surv.edge.csv',row.names = FALSE,quote=FALSE)
write.csv(b, 'reports/thesis/cyto/actf.direct.2.surv.node.csv',row.names = FALSE,quote = FALSE)
##############################

