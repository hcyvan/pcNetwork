source('./lib/helpers.R')
source('./R/lib.R')

diff.gene <- pf.get.diff(c('lncRNA','pcg'))

genes.zfpkm <- pf.filter.zfpkm(diff.gene$GeneID)
cor.pairs <- multicor(genes.zfpkm, rds = './data/cor.pearson.lncRNA.pcg.5343.zfpkm.pairs.rds')



