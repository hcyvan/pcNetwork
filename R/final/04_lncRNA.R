source('./R/lib.R')

lncRNA<-pf.get.diff('lncRNA')
diff<-pf.get.diff()

lncRNA.zfpkm<-pf.filter.zfpkm(lncRNA$GeneID)
zfpkm<-pf.filter.zfpkm(diff$GeneID)

cor.lncRNA.diff <- pf.multicor(t(lncRNA.zfpkm), t(zfpkm), rds = './data/cor.pearson.lncRNA.2.diff5343.zfpkm.pairs.rds')
a<-filter(cor.lncRNA.diff,FDR<0.01,abs(r)>=0.9)
length(unique(as.vector(a$v1)))
