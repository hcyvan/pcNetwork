source('./lib/helpers.R')
source('./R/lib.R')

diff.gene <- pf.get.diff(c('lncRNA','pcg'))

genes.zfpkm <- pf.filter.zfpkm(diff.gene$GeneID)
cor.pairs <- multicor(genes.zfpkm, rds = './data/cor.pearson.lncRNA.pcg.5343.zfpkm.pairs.rds')
cor.pairs2 <- multicor(genes.zfpkm, rds = './data/cor.pearson.lncRNA.pcg.5343.zfpkm.pairs.v2.rds')

cor.pairs.filter<-filter(cor.pairs,FDR<0.01,abs(r)>=0.5)
cor.pairs.type<-mutate(cor.pairs.filter,type1=pf.ensembl2biotype(cor.pairs.filter$v1),type2=pf.ensembl2biotype(cor.pairs.filter$v2), symbol1=pf.ensembl2symbol(cor.pairs.filter$v1),symbol2=pf.ensembl2symbol(cor.pairs.filter$v2))


cor.pairs.1 <- filter(cor.pairs.type, type1%in%pv.lncRNA, type2%in%pv.pcg)
cor.pairs.2 <- filter(cor.pairs.type, type1%in%pv.pcg, type2%in%pv.lncRNA)

cor.pairs.all<-rbind(cor.pairs.1,cor.pairs.2)

write.csv(cor.pairs.all,'./reports/correlation.csv',row.names = FALSE)
