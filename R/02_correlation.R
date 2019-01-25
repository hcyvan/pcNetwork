source('./lib/helpers.R')
source('./R/lib.R')

diff.gene <- pf.get.diff(c('lncRNA','pcg'))
genes.fpkm <- pf.filter.fpkm(diff.gene$GeneID)

cor.pairs <- multicor(genes.fpkm, rds = './cache/cor.pearson.pairs.rds')
