source('./lib/helpers.R')
source('./R/lib.R')

diff.gene <- pf.get.diff(c('lncRNA','pcg'))

genes.fpkm <- pf.filter.fpkm(diff.gene$GeneID)
cor.pairs <- multicor(genes.fpkm, rds = './cache/cor.pearson.pairs.rds')


genes.zfpkm <- pf.filter.zfpkm(diff.gene$GeneID)
cor.pairs.2 <- multicor(genes.zfpkm, rds = './cache/cor.pearson.pairs.2.rds')


